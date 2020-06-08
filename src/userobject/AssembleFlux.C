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
  params.addParam<std::string>("dc_boundaries", "-1", "Dirichlet Boundary ID");
  params.addParam<std::string>("dc_variables" , "-1", "Variable to which given BC_id applies");
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
_qrule(_assembly.qRule()),
_dc_var(getParam<std::string>("dc_variables"))

{
// std::cout<<"ciao IN"<<std::endl;
  if (_hasMeshModifier)
  {
    _meshModifierName=getParam<std::string>("fractureMeshModifier");
    //exit(1);
  }
    
  std::vector<std::string> tmp = split_string(parameters.get<std::string>("dc_boundaries"), ' ');
  for(auto str_tmp=tmp.begin(); str_tmp != tmp.end(); str_tmp++)
  {
        _dc_boundary_id.push_back(atoi(str_tmp->c_str()));
  }
    
// std::cout<<"ciao OUT"<<std::endl;
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
  DenseMatrix<Number> ke_1;
  DenseMatrix<Number> ke_2;
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

          if (elem->subdomain_id()==1)
            ke_1(i,j) += -1.0*JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );
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
     //_diri_flux.print_matlab("diri_flux.m");
     //_neum_flux.print_matlab("neum_flux.m");
     _diri_flux.add(1.0,_neum_flux);
     _tot_flux.add(-1.0,_diri_flux);
     auto _f_tot = _flux_1.sum();
     std::cout<<"f_tot= "<<_f_tot<<std::endl;
    
     _stiffness_matrix_1.add(1.0,*_stiffness_matrix_t);
    

     _stiffness_matrix_1.vector_mult(_flux_1,*_nl.currentSolution());
     _flux_1.add(1.0,_neum_flux);
     _flux_1.add(-1.0,_diri_flux);



     _stiffness_matrix_2.vector_mult(_flux_2,*_nl.currentSolution());

     _flux_1.print_matlab("f_1.m");
     _flux_2.print_matlab("f_2.m");
     auto _f1 = _flux_1.sum();
     auto _f2 = _flux_2.sum();

     std::cout<<"f_1= "<<_f1<<std::endl;
     std::cout<<"f_2= "<<_f2<<std::endl;
    
     AssembleFlux::determine_dc_bnd_var_id(AssembleFlux::split_string(_dc_var, ' '));
    
     ConstBndNodeRange & bnd_nodes = *_fe_problem.mesh().getBoundaryNodeRange();
    
     unsigned int i = 0;

     PetscVector<Number> _bc_vec(_fe_problem.es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
     _bc_vec.zero();
 
     for(auto boundary = _dc_boundary_id.begin(); boundary != _dc_boundary_id.end(); ++boundary, i++)
     {
        // iterate just over boundary nodes
        for (const auto & bnode : bnd_nodes)
        {
            libMesh::Node * current_node = bnode->_node;
            
            // check if node is in active boundary list
            if (_fe_problem.mesh().isBoundaryNode(current_node->id(), *boundary))
            {
                for (auto v = 0; v < _fe_problem.getNonlinearSystemBase().nVariables(); v++)
                {
                    const Variable & var = _nl.system().variable(v);
                    
                    unsigned int var_num = var.number();
                   
                    if (current_node->n_dofs(_fe_problem.getNonlinearSystemBase().number(), var_num) > 0)
                    {
                        
                        if(std::find(_dc_variables_id[i].begin(), _dc_variables_id[i].end(), var_num) != _dc_variables_id[i].end())
                        {
                            _bc_vec.set(current_node->dof_number(_fe_problem.getNonlinearSystemBase().number(), var_num, 0), 1.0);
                            
                        }
                    }
                }
            }
        }
     }
    
     _bc_vec.close();
    

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





    

    



void
AssembleFlux::determine_dc_bnd_var_id(const std::vector<std::string> & BC_var){
    // automatic fill-in
    NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();
    
    std::vector<int> vec(_nl.nVariables());
    
    std::iota(vec.begin(), vec.end(), 0);
    
    unsigned int i;
    
    auto str_tmp = BC_var.begin();
    
    PetscFunctionBegin;
    // going over all BC_ids
    for(i = 0; str_tmp != BC_var.end(); i++, str_tmp++)
    {
        std::vector<std::string> tmp = AssembleFlux::split_string(*str_tmp, '-');
        
        // check if variable assigned in the input file exists for given simulation
        bool var_flg = 1;
        for(auto t = tmp.begin(); t != tmp.end(); ++t)
        {
            if(atoi(t->c_str()) >= _nl.nVariables())
                var_flg = 0;
        }
        
        // in case u havent put anything into input file, or u put too much
        if(*str_tmp == "-1" || var_flg == 0)
        {
            //std::cout<<"no_si"<<_nl.nVariables()<<std::endl;
            _dc_variables_id.push_back(vec);
        }
        else
        {
            unsigned int j;
            std::vector<int > one_BC_id;
            auto str_in = tmp.begin();
            for(j = 0; str_in != tmp.end(); j ++, str_in++)
            {
                one_BC_id.push_back(atoi(str_in->c_str()));
            }
            _dc_variables_id.push_back(one_BC_id);
        }
    }
    
    // check if u have same number of BC_ids in both parameters
    if(_dc_variables_id.size() != _dc_boundary_id.size())
    {
        _dc_variables_id.clear();
        for(auto i = 0; i != _dc_boundary_id.size(); i++)
        {
            _dc_variables_id.push_back(vec);
        }
    }
    
    // print out what is considered for zero-ing
    std::cout<<" ------ BC CONDITIONS Begin ------ \n";
    unsigned int t = 0;
    //std::cout<<"_dc_variables_id.begin()"<<_dc_variables_id.size()<<std::endl;
    for(auto i = _dc_variables_id.begin(); i != _dc_variables_id.end();  t++, i++)
    {
        std::cout<<"\n BC_id:  "<< _dc_boundary_id[t] << "   var_ids:  ";
        std::for_each(i->begin(), i->end(), [](int i){ std::cout << i << "  " ; });
    }
    std::cout<<" ------ BC CONDITIONS End ------ \n";
}





std::vector<std::string>
AssembleFlux::split_string(const std::string & s, char delim)
{
    
    std::vector<std::string> v;
    
    if (s.length() == 0)
        std::cerr << "Got an empty string. Split_string(...) is confused. \n";
    
    auto i = 0;
    auto pos = s.find(delim);
    while (pos != std::string::npos)
    {
        v.push_back(s.substr(i, pos - i));
        i = ++pos;
        pos = s.find(delim, pos);
        
        if (pos == std::string::npos)
            v.push_back(s.substr(i, s.length()));
    }
    
    if (v.size() == 0) // if only one word is in the string
        v.push_back(s);
    
    return v;
}

