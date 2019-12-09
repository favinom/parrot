//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolveDiffusion.h"
#include <iostream>
#include <string>
#include "FEProblem.h"
#include "MooseVariableFEBase.h"
#include "NonlinearSystemBase.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh_base.h"
//#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"
//#include "libmesh/nonlinear_implicit_system.h"
//#include "libmesh/transient_system.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exodusII_io.h"


#include "libmesh/quadrature_gauss.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/parsed_function.h"

using namespace std;

registerMooseObject("parrotApp", SolveDiffusion);


template <>
InputParameters
validParams<SolveDiffusion>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addRequiredParam<std::vector<int>>("block_id","block_id");
  params.addRequiredParam<std::vector<Real>>("value_p","value_p");
  params.addRequiredParam<std::vector<boundary_id_type>>("boundary_D_bc","boundary_D_bc");
  params.addRequiredParam<std::vector<boundary_id_type>>("boundary_N_bc","boundary_N_bc");
  params.addRequiredParam<std::vector<Real>>("value_N_bc", "The value of Neumann");
  params.addRequiredParam<std::vector<Real>>("value_D_bc", "The value of Dirichlet");
  params.addRequiredParam<std::vector<AuxVariableName>>("aux_variable", "The auxiliary variable to store the transferred values in.");
  params.addParam<std::string>("output_file", "the file name of the output");
  params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");

  return params;
}

SolveDiffusion::SolveDiffusion(const InputParameters & parameters) :
GeneralUserObject(parameters),
_aux_var_names(getParam<std::vector<AuxVariableName>>("aux_variable")),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
_boundary_D_ids(getParam<std::vector<boundary_id_type>>("boundary_D_bc")),
_boundary_N_ids(getParam<std::vector<boundary_id_type>>("boundary_N_bc")),
_value_N_bc(getParam<std::vector<Real>>("value_N_bc")),
_value_D_bc(getParam<std::vector<Real>>("value_D_bc")),
_has_output_file( isParamValid("output_file") ),
_hasMeshModifier( isParamValid("fractureMeshModifier") )
{

  if (_has_output_file)
    _output_filename=getParam<std::string>("output_file");

  if (_hasMeshModifier)
    _meshModifierName=getParam<std::string>("fractureMeshModifier");

  _sys_name="Diffusion";
  _var_name="pressure";

  if (_aux_var_names.size() == 1)
    _aux_var_name = _aux_var_names.at(0);
  else
    paramError("variable", "You need to specify one and only one variable");

}



void SolveDiffusion::initialize()
{
  _console<<"BEGIN initialize\n";

  EquationSystems & mooseEquationSystems=_fe_problem.es();
   // get_mesh gives a reference to a meshbase
  MeshBase & mooseMesh=mooseEquationSystems.get_mesh();
  EquationSystems equation_systems ( mooseMesh );
  LinearImplicitSystem & system = equation_systems.add_system<LinearImplicitSystem> (_sys_name.c_str());
  _p_var = system.add_variable (_var_name.c_str(), FIRST);
  equation_systems.init();
  equation_systems.print_info();

  solve(equation_systems);
  mooseEquationSystems.reinit();
  set_solution(equation_systems);

  equation_systems.delete_system(_sys_name.c_str());
  equation_systems.clear();

  _console<<"END initialize\n";
}


int SolveDiffusion::solve(EquationSystems & es)
{
  _console<<"BEGIN solve\n";

  LinearImplicitSystem & _system = es.get_system<LinearImplicitSystem> (_sys_name);

  AssembleDiffusionOP(es, _sys_name.c_str() );

  NumericVector<Number> & sol_NV=*_system.solution;
  NumericVector<Number> & rhs_NV=*_system.rhs;
  SparseMatrix<Number>  & mat_SM=*_system.matrix;

  PetscVector<Number> & sol_PV= dynamic_cast<PetscVector<Number> &>(sol_NV);
  PetscVector<Number> & rhs_PV= dynamic_cast<PetscVector<Number> &>(rhs_NV);
  PetscMatrix<Number> & mat_PM= dynamic_cast<PetscMatrix<Number> &>(mat_SM);

  PetscErrorCode ierr;
  PC _diff_problem;

  ierr = PCCreate(PETSC_COMM_WORLD, &_diff_problem);
  CHKERRQ(ierr);
  ierr = PCSetType(_diff_problem,PCLU);
  CHKERRQ(ierr);
  ierr = PCSetOperators(_diff_problem, mat_PM.mat(),mat_PM.mat());
  CHKERRQ(ierr);  
  ierr = PCFactorSetMatSolverPackage(_diff_problem,MATSOLVERMUMPS);
  CHKERRQ(ierr);
  ierr = PCApply(_diff_problem,rhs_PV.vec(),sol_PV.vec()); CHKERRQ(ierr);
  CHKERRQ(ierr);
  PCDestroy(&_diff_problem);
  //solution->print_matlab();

  if(_has_output_file)
    ExodusII_IO (es.get_mesh()).write_equation_systems(_output_filename.c_str(), es);

  _console<<"END solve\n";

  return 0 ;
}

void SolveDiffusion::AssembleDiffusionOP(EquationSystems & _es, std::string const & system_name)
{
  _console << "BEGIN Assemble_Diffusion"  << std::endl;

  MeshModifier       const * _myMeshModifier_ptr;
  FractureUserObject const * _fractureUserObject_ptr;

  if (_hasMeshModifier)
  {
  MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
  _myMeshModifier_ptr=&_myMeshModifier;
  FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
  _fractureUserObject_ptr=&_fractureUserObject;
}
    // Get a constant reference to the mesh object.
  const MeshBase & mesh = _es.get_mesh();
    // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();
    // Get a reference to our system.
  LinearImplicitSystem & _system = _es.get_system<LinearImplicitSystem>(system_name);

  SparseMatrix<Number> & matrix_K = *_system.matrix;
  matrix_K.zero();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
  FEType fe_type = _system.get_dof_map().variable_type(0);

  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

  QGauss qrule (dim, TENTH );
    //QGauss qrule (dim, fe_type.default_quadrature_order() );
    // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));

    // Boundary integration requires another quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    // In 1D, the Clough and Gauss quadrature rules are identical.
  std::unique_ptr<QBase> qface(fe_type.default_quadrature_rule(dim-1));

    // Tell the finite element object to use our
    // quadrature rule.
  fe_face->attach_quadrature_rule (qface.get());

    // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
  const std::vector<Point> & q_points = fe->get_xyz();
  const DofMap & dof_map = _system.get_dof_map();



  std::vector<dof_id_type> dof_indices;

  DenseMatrix<Number> ke;
  DenseVector<Number> re;

  MeshBase::const_element_iterator           el = mesh.active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem * elem = *el;

    fe->reinit (elem);
        //std::cout<<qrule.n_points()<<std::endl;

    dof_map.dof_indices(elem, dof_indices);

    const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());

        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
    libmesh_assert_equal_to (n_dofs, phi.size());

        //std::cout<<"n_dofs "<<n_dofs<<std::endl;

    ke.resize (n_dofs , n_dofs);
    ke.zero();

    for (unsigned int i=0; i<phi.size(); i++){

      for (unsigned int j=0; j<phi.size(); j++){

        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

          Real permeabiltiy=ComputeMaterialProprties(elem);

          if(_hasMeshModifier)
          {
            if ( _fractureUserObject_ptr[0].isInside(q_points[qp]) )
            {
             permeabiltiy=_vector_value.at(_vector_value.size()-1);
            }
          }


          ke(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeabiltiy * dphi[i][qp] ) );
        }
      }

    }

    re.resize(n_dofs);
    re.zero();

    
    for (auto side : elem->side_index_range())
    {
      if (elem->neighbor_ptr(side) == nullptr)
      {
        const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
        const std::vector<Real> & JxW_face = fe_face->get_JxW();

        fe_face->reinit(elem, side);

        for(int k =0; k < _boundary_N_ids.size(); k++){

                if (mesh.get_boundary_info().has_boundary_id (elem, side, _boundary_N_ids[k])) // Apply a traction on the right side
                {
                 for (unsigned int qp=0; qp<qface->n_points(); qp++)
                   for (unsigned int i=0; i != n_dofs; i++)
                    re(i) += JxW_face[qp] * -1.0 * _value_N_bc.at(k) * phi_face[i][qp];
                }
              }
            }
          }

          const Real penalty = 1.e15;

          for (auto s : elem->side_index_range())
          {

            if (elem->neighbor_ptr(s) == nullptr)
            {
              fe_face->reinit(elem, s);

              const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
              //const std::vector<Real> & JxW_face = fe_face->get_JxW();

              std::vector<Real> whichOnBoundary(phi_face.size(),0.0);

              for (int i=0; i<phi_face.size(); ++i)
              {
                for (int qp=0; qp<phi_face.at(i).size(); ++qp)
                  whichOnBoundary.at(i)+=phi_face.at(i).at(qp);
              }

              Real boundaryValue=0.0;

              //std::cout<<"BID" <<mesh.get_boundary_info().boundary_id(elem, s)<<std::endl;

              for(int k =0; k < _boundary_D_ids.size(); k++)
              {

                if (mesh.get_boundary_info().has_boundary_id (elem, s, _boundary_D_ids[k]))
                {
                  //std::cout<<"_boundary_D_ids[k]="<<_boundary_D_ids[k]<<" _value_D_bc.at(k)"<<_value_D_bc.at(k)<<std::endl;
                 boundaryValue=_value_D_bc.at(k);


                   //_console << "Assemble_Diffusion:: Add penalty BC on Id"  << _boundary_D_ids[k] <<std::endl;


                 for (int k=0; k<whichOnBoundary.size(); ++k)
                 {
                  if (whichOnBoundary.at(k)>0.5)
                  {
                   re(k) += penalty * boundaryValue;
                   ke(k,k) += penalty;
                 }
               }
             }
           }
         }
       }

       dof_map.constrain_element_matrix_and_vector (ke, re, dof_indices);

       _system.matrix->add_matrix (ke, dof_indices);
       _system.rhs->add_vector(re, dof_indices);
     }



     _system.matrix->close();
     _system.rhs->close();
     _console << "END Assemble_Diffusion"  << std::endl;    
   }


   Real
   SolveDiffusion::ComputeMaterialProprties(const Elem *elem)
   {
    Real permeability=0.0;
    for(int ll=0; ll<_vector_p.size(); ll++)
    {
      if (elem->subdomain_id()==_vector_p[ll])
      {
        permeability = _vector_value[ll];
      }
    }
    return permeability;
  }




  void SolveDiffusion::set_solution(EquationSystems & _es)
  {
    MeshBase & _mesh = _fe_problem.es().get_mesh();


    
    MooseVariableFEBase  & aux_var = _fe_problem.getVariable(0, _aux_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
    
    // solution of the original system
    System & aux_sys = aux_var.sys().system();



    NumericVector<Number> * aux_solution = aux_sys.solution.get();

    //NumericVector<Number> & from_solution = *ls.solution;
    

    LinearImplicitSystem & ls = _es.get_system<LinearImplicitSystem>("Diffusion");


    { // loop through our local elements and set the solution from the projection

      for (const auto & node : _mesh.local_node_ptr_range())

      {
        for (unsigned int comp = 0; comp < node->n_comp(aux_sys.number(), aux_var.number()); comp++)

        {

            //std::cout<<"uno"<<std::endl;

          const dof_id_type proj_index = node->dof_number(ls.number(), _p_var, comp);

            //std::cout<<"due"<<std::endl;

          const dof_id_type to_index = node->dof_number(aux_sys.number(), aux_var.number(), comp);

            //std::cout<<"tre"<<std::endl;

             //

          aux_solution->set(to_index, (*ls.solution)(proj_index));
        }

      }
    }

    aux_solution->close();
    aux_sys.update();

  }
