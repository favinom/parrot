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
#include "libmesh/quadrature_gauss.h"
#include "FEProblem.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/parsed_function.h"
#include <iostream> 
#include <string> 
using namespace std;

registerMooseObject("parrotApp", SolveDiffusion);

void
assembly_diffusion(EquationSystems & es, const std::string & system_name)
{
    std::cout<<"CIAO\n";
    SolveDiffusion * pippo = es.parameters.get< SolveDiffusion*>("pippo");
    pippo->AssembleDiffusionOP(es, system_name);
}


template <>
InputParameters
validParams<SolveDiffusion>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addRequiredParam<std::vector<int>>("block_id",
                                                     "The name of the nodeset to create");
  params.addRequiredParam<std::vector<Real>>("value_p",
                                                     "The name of the nodeset to create");

  params.addRequiredParam<std::vector<boundary_id_type>>("boundary_D_bc",
                                                     "The name of the boundary to create");

  params.addRequiredParam<std::vector<boundary_id_type>>("boundary_N_bc",
                                                     "The name of the boundary to create");

  params.addRequiredParam<std::vector<std::string>>("function_D_bc",
                                                     "The value of the dirichlet_bc");

  params.addRequiredParam<std::vector<Real>>("value_N_bc", "The vlue of Neumann");


  return params;
}

SolveDiffusion::SolveDiffusion(const InputParameters & parameters) :
GeneralUserObject(parameters),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
_boundary_D_ids(getParam<std::vector<boundary_id_type>>("boundary_D_bc")),
_boundary_N_ids(getParam<std::vector<boundary_id_type>>("boundary_N_bc")),
_string_expr(getParam<std::vector<std::string>>("function_D_bc")),
_value_N_bc(getParam<std::vector<Real>>("value_N_bc"))

{
    
   
}

 void SolveDiffusion::execute(){};


 void SolveDiffusion::initialize()
{
    std::cout<<"initialize\n";

   // Create an equation systems object.
   EquationSystems equation_systems (_fe_problem.es().get_mesh());

   LinearImplicitSystem & _system = equation_systems.add_system<LinearImplicitSystem> ("Diffusion");

   equation_systems.parameters.set<SolveDiffusion *>("pippo") = this;

   _p_var = _system.add_variable ("pressure", FIRST);

   _system.attach_assemble_function(assembly_diffusion);

   add_bc(_system);

   equation_systems.init();

   equation_systems.print_info();

   solve(equation_systems, _system);

   set_solution(equation_systems);


   //    _p_var = _fe_problem.es().add_system<LinearImplicitSystem>("Diffusion").add_variable("pressure",FIRST),
    
   //    _fe_problem.es().reinit();
 
    
}

 void SolveDiffusion::finalize()
{
    std::cout<<"finalize\n";
}

 void SolveDiffusion::threadJoin(const UserObject & uo)
{
	
};



void SolveDiffusion::solve(EquationSystems & _es, LinearImplicitSystem & _system){


	   std::cout<<"solve\n";

	   MeshBase & _mesh_m = _fe_problem.es().get_mesh();

       //LinearImplicitSystem & _system = equation_systems.get_system<LinearImplicitSystem> ("Diffusion");
       
       libMesh::Parallel::Communicator const & comm_in=_mesh_m.comm();

       _system.solve();
    
       NumericVector<Number> * solution = _system.solution.get();

       //solution->print_matlab();

       ExodusII_IO (_es.get_mesh()).write_equation_systems("matrix_c.e", _es);
       
       NumericVector<Number> & sol_ref = solution[0];

}

void SolveDiffusion::AssembleDiffusionOP(EquationSystems & _es, const std::string & system_name){

	_console << "Assemble_Diffusion::Begin "  << std::endl;

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _es.get_mesh();

    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to our system.
    LinearImplicitSystem & _system = _es.get_system<LinearImplicitSystem>(system_name);

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);

    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

    SparseMatrix<Number> & matrix_K = *_system.matrix;
 
    DenseMatrix<Number> ke;

    DenseVector<Number> re;

    QGauss qrule (dim, fe_type.default_quadrature_order());

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

    const DofMap & dof_map = _system.get_dof_map();

    std::vector<dof_id_type> dof_indices;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // first we need to manually zero the matrix
    matrix_K.zero();

    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;

        Elem * ele = *el;

        fe->reinit (elem);

        dof_map.dof_indices(elem, dof_indices);

        const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());

        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        libmesh_assert_equal_to (n_dofs, phi.size());

        //std::cout<<"n_dofs "<<n_dofs<<std::endl;

        ke.resize (n_dofs , n_dofs);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                    ke(i,j) += ComputeMaterialProprties(elem) * JxW[qp] * dphi[i][qp] * dphi[j][qp];
                }
            }

        }

        re.resize(n_dofs);

    {
        for (auto side : elem->side_index_range())
        {
          if (elem->neighbor_ptr(side) == nullptr)
            {
              const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              fe_face->reinit(elem, side);

              for(int k =0; k < _boundary_N_ids.size(); k++){

                std::cout<<_boundary_N_ids[k] <<" "<<_value_N_bc.at(k) << std::endl;

                if (mesh.get_boundary_info().has_boundary_id (elem, side, _boundary_N_ids[k])) // Apply a traction on the right side
                  {
                     for (unsigned int qp=0; qp<qface->n_points(); qp++)
                     for (unsigned int i=0; i != n_dofs; i++)
                     		re(i) += JxW_face[qp] * -1.0 * _value_N_bc.at(k)  * phi_face[i][qp];
                  }
               }
             }
          }
     }




        dof_map.heterogenously_constrain_element_matrix_and_vector (ke, re, dof_indices);

        _system.matrix->add_matrix (ke, dof_indices);

        _system.rhs->add_vector(re, dof_indices);

        //exit(1);



    }



    _system.matrix->close();

    _system.matrix->print_matlab("Diff.m");

    _console << "Assemble_Diffusion::end "  << std::endl;

    
}


Real
SolveDiffusion::ComputeMaterialProprties(const Elem *elem){
    
   // _console << "_vector_p.size()"  << _vector_p.size() <<std::endl;

Real permeability=0.0;

    for(int ll=0; ll<_vector_p.size(); ll++){
        if (elem->subdomain_id()==_vector_p[ll]) {

            permeability = _vector_value[ll]; 
        }
    }
  
    return permeability;   
}

void SolveDiffusion::add_bc(LinearImplicitSystem & _system){

	_console << "SolveDiffusion::add_bc"  << std::endl;

  std::vector<ParsedFunction<Number> > pf;

  pf.clear();

  DofMap & dof_map = _system.get_dof_map();

  // ParsedFunction<Number> pf_d2(_string_expr);

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor

  std::vector<unsigned int> variables(1);

  variables[0] = _p_var; 

  for(int k =0; k < _boundary_D_ids.size(); k++){

  	    std::set<boundary_id_type> boundary_dr;

  	    boundary_dr.clear();

        boundary_dr.insert(_boundary_D_ids[k]);

        ParsedFunction<Number> temp(_string_expr.at(k));

        pf.push_back(temp);

        std::cout<<_boundary_D_ids[k] <<" "<<_string_expr.at(k) << std::endl;

        dof_map.add_dirichlet_boundary(DirichletBoundary(boundary_dr,
                                                         variables,
                                                         pf.at(k), LOCAL_VARIABLE_ORDER));

  }



  //  DirichletBoundary dirichlet_bc(boundary_dr, variables, pf,
                                 //LOCAL_VARIABLE_ORDER);


  //LinearImplicitSystem & _system = _fe_problem.es().get_system<LinearImplicitSystem>("Diffusion");

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  //_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

}



void SolveDiffusion::set_solution(EquationSystems & _es)
{
    // copy projected solution into target es
    MeshBase & _mesh = _fe_problem.es().get_mesh();
    
    MooseVariable & _var = _fe_problem.getStandardVariable(0, "pressure");
    
    // solution of the original system
    System & _sys = _var.sys().system();

    NumericVector<Number> * _solution = _sys.solution.get();
    
    LinearImplicitSystem & ls = _es.get_system<LinearImplicitSystem>("Diffusion");
    { // loop through our local elements and set the solution from the projection
        
        for (const auto & node : _mesh.local_node_ptr_range())

        {
        	for (unsigned int comp = 0; comp < node->n_comp(_sys.number(), _var.number()); comp++)

        	{
        		const dof_id_type proj_index = node->dof_number(ls.number(), _p_var, comp);

        		const dof_id_type to_index = node->dof_number(_sys.number(), _var.number(), comp);

        		_solution->set(to_index, (*ls.solution)(proj_index));
        	}

        }
    }
    
    
    _solution->close();

    _sys.update();
}
