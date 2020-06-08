//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssembleFlux.h"
#include <iostream>
#include <string>
#include "FEProblem.h"
#include "MooseVariableFEBase.h"
#include "Assembly.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"

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
#include "libmesh/nemesis_io.h"


#include "libmesh/quadrature_gauss.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/quadrature_grid.h"
#include "libmesh/quadrature.h"

using namespace std;

registerMooseObject("parrotApp", AssembleFlux);


template <>
InputParameters
validParams<AssembleFlux>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addRequiredParam<std::vector<int>>("block_id","block_id");
  params.addRequiredParam<std::vector<Real>>("value_p","value_p");
  params.addRequiredParam<std::vector<boundary_id_type>>("boundary_D_bc","boundary_D_bc");
  params.addRequiredParam<std::vector<boundary_id_type>>("boundary_N_bc","boundary_N_bc");
  params.addRequiredParam<std::vector<Real>>("value_N_bc", "The value of Neumann");
  params.addRequiredParam<std::vector<Real>>("value_D_bc", "The value of Dirichlet");
  params.addRequiredParam<AuxVariableName>("sol_variable", "The auxiliary variable to store the transferred values in.");
  //params.addParam<std::string>("output_file", "the file name of the output");
  params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  params.addRequiredParam<bool>("conservative","use a conservative scheme?");
//  params.addRequiredParam<int>("solver_type","solver_type");

  return params;
}


AssembleFlux::AssembleFlux(const InputParameters & parameters) :
GeneralUserObject(parameters),
//_pp_comm(_fe_problem.es().get_mesh().comm()),
//_stiffness_matrix_1(_pp_comm),
//_stiffness_matrix_2(_pp_comm),
//_stiffness_matrix_t(_pp_comm),
_sol_var_name(getParam<AuxVariableName>("sol_variable")),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
_boundary_D_ids(getParam<std::vector<boundary_id_type>>("boundary_D_bc")),
_boundary_N_ids(getParam<std::vector<boundary_id_type>>("boundary_N_bc")),
_value_N_bc(getParam<std::vector<Real>>("value_N_bc")),
_value_D_bc(getParam<std::vector<Real>>("value_D_bc")),
_hasMeshModifier( isParamValid("fractureMeshModifier") ),
_conservativeScheme( getParam<bool>("conservative") ),
_qrule(_assembly.qRule())
{
 std::cout<<"ciao IN"<<std::endl;
  if (_hasMeshModifier)
  {
    _meshModifierName=getParam<std::string>("fractureMeshModifier");
    //exit(1);
  }
 std::cout<<"ciao OUT"<<std::endl;
 }



void
AssembleFlux::initialize()
{
  _console<<"BEGIN initialize\n";
    
  DofMap   const & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();

//  _stiffness_matrix_1.attach_dof_map(dof_map);
//  _stiffness_matrix_1.zero();
//
//  _stiffness_matrix_2.attach_dof_map(dof_map);
//  _stiffness_matrix_2.zero();
//
//  _stiffness_matrix_t.attach_dof_map(dof_map);
//  _stiffness_matrix_t.zero();
    
  solve();
    
  _console<<"END initialize\n";
}



int
AssembleFlux::solve()
{
  _console<<"BEGIN Assembly\n";

  ComputeFlux();
  bool assemble=true;

  _console<<"END Assembly\n";
    return 1;
 }

void
AssembleFlux::ComputeFlux()
{
  _console << "BEGIN ComputeFlux"  << std::endl;

  MeshModifier       const * _myMeshModifier_ptr;
  FractureUserObject const * _fractureUserObject_ptr;
    
  DofMap   const & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();
    
  PetscMatrix<Number> _stiffness_matrix_1(_fe_problem.es().get_mesh().comm());
  PetscMatrix<Number> _stiffness_matrix_2(_fe_problem.es().get_mesh().comm());
  //PetscMatrix<Number> _stiffness_matrix_t(_fe_problem.es().get_mesh().comm());
    
  _stiffness_matrix_1.attach_dof_map(dof_map);
  _stiffness_matrix_1.init();
  //_stiffness_matrix_1.clear();

  _stiffness_matrix_2.attach_dof_map(dof_map);
  _stiffness_matrix_2.init();
  //_stiffness_matrix_2.clear();

  //_stiffness_matrix_t.attach_dof_map(dof_map);
  //_stiffness_matrix_t.init();
  //_stiffness_matrix_t.clear();


  if (_hasMeshModifier)
  {
    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    _myMeshModifier_ptr=&_myMeshModifier;
    FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
    _fractureUserObject_ptr=&_fractureUserObject;
  }

  const MeshBase & mesh = _fe_problem.es().get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  PetscVector<Number> _neum_flux(_fe_problem.es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
  _neum_flux.zero();

  PetscVector<Number> _diri_flux(_fe_problem.es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
  _diri_flux.zero();
    
  PetscVector<Number> _tot_flux(_fe_problem.es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
  _tot_flux.zero();

  PetscVector<Number> _flux_1(_fe_problem.es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
  _flux_1.zero();

  PetscVector<Number> _flux_2(_fe_problem.es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
  _flux_2.zero();
    
  

  TransientNonlinearImplicitSystem const & _system = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
  FEType fe_type = _system.get_dof_map().variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  std::unique_ptr<QBase> qrule( QBase::build (_qrule->type(),dim,_qrule->get_order()));
  fe->attach_quadrature_rule (qrule.get());
  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
  std::unique_ptr<QBase> qface(fe_type.default_quadrature_rule(dim-1));
  fe_face->attach_quadrature_rule (qface.get());


  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
  const std::vector<Point> & q_points = fe->get_xyz();
  //const DofMap & dof_map = _system.get_dof_map();



  std::vector<dof_id_type> dof_indices;

  DenseMatrix<Number> ke_1, ke_2;
  DenseMatrix<Number> ke_t;
  DenseVector<Number> re;

  MeshBase::const_element_iterator           el = mesh.active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();

  std::vector<Number> permeability;
 

  for ( ; el != end_el; ++el)
  {
    const Elem * elem = *el;
    fe->reinit (elem);

    dof_map.dof_indices(elem, dof_indices);
    const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());

    libmesh_assert_equal_to (n_dofs, phi.size());

    ke_1.resize (n_dofs , n_dofs);
    ke_1.zero();

    ke_2.resize (n_dofs , n_dofs);
    ke_2.zero();

    //ke_t.resize (n_dofs , n_dofs);
    //ke_t.zero();

    Real localPermeability=ComputeMaterialProprties(elem);
    permeability.assign( qrule->n_points() ,localPermeability );

    if(_hasMeshModifier)
    {
      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      {
        if ( _fractureUserObject_ptr[0].isInside(q_points[qp]) )
        {
          permeability.at(qp)=_vector_value.at(_vector_value.size()-1);
        }
      }
    }

    for (unsigned int i=0; i<phi.size(); i++)
      for (unsigned int j=0; j<phi.size(); j++)
        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
         //ke_t(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );

          if (elem->subdomain_id()==0)
            ke_1(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );
          if (elem->subdomain_id()==1)
            ke_2(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );
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

       dof_map.constrain_element_matrix_and_vector (ke_1, re, dof_indices);
       dof_map.constrain_element_matrix (ke_2,dof_indices);
       //dof_map.constrain_element_matrix (ke_t,dof_indices);
  
       _stiffness_matrix_1.add_matrix (ke_1, dof_indices);
       _stiffness_matrix_2.add_matrix (ke_2, dof_indices);
       //_stiffness_matrix_t.add_matrix (ke_t, dof_indices);

       _neum_flux.add_vector(re, dof_indices);
     }



     MooseVariable & var = _fe_problem.getStandardVariable(0, _sol_var_name);

     System & sys = var.sys().system();

     auto _sol=sys.current_local_solution.get();

     _stiffness_matrix_1.close();
     _stiffness_matrix_2.close();
    //_stiffness_matrix_t.close();
     _neum_flux.close();
    
     auto &_sys = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
     libMesh::PetscMatrix<libMesh::Number> * _stiffness_matrix_t  = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_sys.matrix);
    
     NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();
    
     _fe_problem.computeJacobianSys(_sys, *_nl.currentSolution(), *_stiffness_matrix_t);


     _stiffness_matrix_t->vector_mult(_diri_flux,*_nl.currentSolution());

     _diri_flux.add(-1.0,_neum_flux);
     _tot_flux.add(-1.0,_diri_flux);
     auto _f_tot = _flux_1.sum();
     std::cout<<"f_tot= "<<_f_tot<<std::endl;

     _stiffness_matrix_1.vector_mult(_flux_1,*_nl.currentSolution());
     _flux_1.add(-1.0,_neum_flux);
     _flux_1.add(-1.0,_diri_flux);
//     _flux_1.add(1.0,_diri_flux);


     _stiffness_matrix_2.vector_mult(_flux_2,*_nl.currentSolution());

     _flux_1.print_matlab("f_1.m");
     _flux_2.print_matlab("f_2.m");
     auto _f1 = _flux_1.sum();
//     auto _f2 = _flux_2.sum();
//
     std::cout<<"f_1= "<<_f1<<std::endl;
//     std::cout<<"f_2= "<<_f2<<std::endl;

     _console << "END ComputeFlux"  << std::endl;
}


Real
AssembleFlux::ComputeMaterialProprties(const Elem *elem)
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





