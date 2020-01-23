/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*    Immersed_Boundary- ICS Mechanical simulation framework    */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/



#include "HybridModel.h"
#include "Parrot_preonly.h"
#include "ksp_parrot_impl.h"
#include "libmesh/petsc_nonlinear_solver.h"


#include <petscksp.h>
#include <petscmat.h>

#include "utopia.hpp"
#include "utopia_Socket.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_NewTransferAssembler.hpp"

#include <queue>
#include <iostream>
#include <string>
#include <algorithm> 

#include "MultiApp.h"
#include "FEProblem.h"
#include "Transient.h"
#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseVariableFE.h"
#include "NonlinearSystemBase.h"
#include "Assembly.h"
#include "TransientMultiApp.h"
#include "FEProblemBase.h"
#include "DisplacedProblem.h"
#include "AddVariableAction.h"
#include "NonlinearSystemBase.h"
#include "MooseVariable.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh_base.h"
#include "libmesh/linear_implicit_system.h"
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
#include "libmesh/quadrature_grid.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/boundary_mesh.h"

#include "utopia_LibMeshBackend.hpp"


using namespace std;

registerMooseObject("parrotApp",HybridModel);


typedef utopia::USparseMatrix SparseMatT;
typedef utopia::UVector VecT;




template <>
InputParameters
validParams<HybridModel>()
{
    InputParameters params = validParams<GeneralUserObject>();
    params.addParam<std::string>("dc_boundaries_m", "-1", "Dirichlet Boundary ID");
    params.addParam<std::string>("dc_variables_m" , "-1", "Variable to which given BC_id applies");
    
    params.addRequiredParam<std::vector<int>>("block_id",
                                                     "The name of the nodeset to create");
    params.addRequiredParam<std::vector<Real>>("value_p",
                                                     "The name of the nodeset to create");
    // // params.addRequiredParam<AuxVariableName>("lagrange_variable",
    // //                                          "The auxiliary variable to store the transferred values in.");
    params.addRequiredParam<VariableName>("matrix_variable",
                                          "The variable to transfer from.");
    params.addRequiredParam<VariableName>("fracture_variable",
                                          "The variable to transfer from.");
    //params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    params.addParam<bool>("pressure","false","put true if you solve pressure system");
    params.addParam<bool>("transport","false","put true if you transport system");
    params.addRequiredParam<MultiAppName>("multi_app", "The MultiApp's name in your input file!");
    
    return params;
}

HybridModel::HybridModel(const InputParameters & parameters) :
GeneralUserObject(parameters),
_dc_var_m(getParam<std::string>("dc_variables_m")),
_f_var_name(getParam<VariableName>("fracture_variable")),
_m_var_name(getParam<VariableName>("matrix_variable")),
// _operator_storage(getUserObject<StoreOperators>("operator_userobject")),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
_pressure(getParam<bool>("pressure")),
_transport(getParam<bool>("transport")),
_multiapp_name(getParam<MultiAppName>("multi_app"))


{
    std::vector<std::string> tmp = split_string(parameters.get<std::string>("dc_boundaries_m"), ' ');
    for(auto str_tmp=tmp.begin(); str_tmp != tmp.end(); str_tmp++)
    {
        _dc_boundary_id_m.push_back(atoi(str_tmp->c_str()));
    }

}





void
HybridModel::initialize()
{
    _console << "Initial Setup of Fracture App " << std::endl;
    
    using namespace utopia;
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    
    FEProblemBase & m_problem = _multi_app.problemBase();
    
    FEProblemBase & f_problem = _multi_app.appProblemBase(0);
    
    MeshBase *_f_mesh  = &f_problem.mesh().getMesh();
    
    MeshBase *_m_mesh  = &m_problem.mesh().getMesh();
    
    MooseVariable & _f_var = f_problem.getStandardVariable(0,_f_var_name);
    
    MooseVariable & _m_var = m_problem.getStandardVariable(0, _m_var_name);


    auto &_m_sys = m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    PetscErrorCode ierr; 

    PetscNonlinearSolver<libMesh::Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<libMesh::Number> *>((m_problem.getNonlinearSystemBase()).nonlinearSolver());


    SNES snes = petsc_solver->snes();    

    ierr = SNESGetKSP(snes,&ksp);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);


    if(!_fe_problem.hasUserObject(_userobject_name_A_tot)){

            std::string class_name = "StoreUtopiaOperators";
            
            auto params = m_problem.getMooseApp().getFactory().getValidParams(class_name);
            
            params.set<bool>("use_displaced_mesh") = false;
            
            params.set<ExecFlagEnum>("execute_on") = "initial";

            _fe_problem.addUserObject("StoreUtopiaOperators", _userobject_name_A_tot, params);

            _fe_problem.addUserObject("StoreUtopiaOperators", _userobject_name_M_f, params);

            _fe_problem.addUserObject("StoreUtopiaOperators", _userobject_name_M_m, params);

    }

    if  (_fe_problem.timeStep()==1){

        T_ = std::make_shared<SparseMatT>();

        NewTransferAssembler transfer_assembler;

        TransferOptions opts;

        opts.from_var_num = _m_var.number();
        opts.to_var_num   = _f_var.number();
        unsigned int n_var = 1 ;
        opts.n_var        = n_var;
        opts.tags         = {};


        transfer_assembler.remove_incomplete_intersections(false);

        transfer_assembler.assemble(
                                    make_ref(*_m_mesh),
                                    make_ref(_m_var.sys().system().get_dof_map()),
                                    make_ref(*_f_mesh),
                                    make_ref(_f_var.sys().system().get_dof_map()), 
                                    opts);

        auto op = transfer_assembler.build_operator();

        T_ = op->matrix();

        set_zero_at_constraint_rows(_m_var.sys().system().get_dof_map(), *T_);



        if (_transport)
        {

            // assemble_mass_matrix(_porosity_m, m_problem,  mass_m, mass_lumped_m);

            // assemble_mass_matrix(_porosity_f, f_problem, mass_f, mass_lumped_f);

            assemble_poro_mass_matrix(m_problem, _userobject_name_M_m);

            assemble_poro_mass_matrix(f_problem, _userobject_name_M_f);
        }
    }


}



bool
HybridModel::solve(){

    bool ok=false;

    if(_pressure) solve_pressure_monolithic(ok);

    if(_transport) solve_transport_monolithic(ok);

    return ok;
    
}


void
HybridModel::execute()
{

    bool done = solve(); 
        
    if (done){
        
        CopyMatrixSolution(sol_m);
        
        CopyFractureSolution(sol_f);
    }


   
}


int HybridModel::solve_pressure_monolithic(bool & ok)
{
    _console << "Solve_pressure_monolithic()"  << std::endl;
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    

    
    FEProblemBase & _f_problem = _multi_app.appProblemBase(0);
    
    MooseVariable & _f_var   = _f_problem.getStandardVariable(0, _f_var_name);
    
    auto &_f_sys = _f_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_f = _f_problem.getNonlinearSystemBase();
    
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_f = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_f_sys.matrix);
    
    _f_problem.computeResidualSys(_f_sys, *_nl_f.currentSolution(), _nl_f.RHS());
    
    _f_problem.computeJacobianSys(_f_sys, *_nl_f.currentSolution(), *petsc_mat_f);

    
    FEProblemBase & _m_problem = _multi_app.problemBase();

    MooseVariable & _m_var = _m_problem.getStandardVariable(0, _m_var_name);

    auto &_m_sys = _m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_m = _m_problem.getNonlinearSystemBase();

    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_m_sys.matrix);
    
    _m_problem.computeResidualSys(_m_sys, *_nl_m.currentSolution(), _nl_m.RHS());
    
    _m_problem.computeJacobianSys(_m_sys, *_nl_m.currentSolution(), *petsc_mat_m);

    utopia::UVector rhs_f;

    utopia::USparseMatrix A_f;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_f.RHS()), rhs_f);

    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_f).mat(), A_f);
    
    utopia::UVector rhs_m;

    utopia::USparseMatrix A_m;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_m.RHS()), rhs_m);

    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_m).mat(), A_m);

    _console << "Solve_monolithic():: START SOLVING"  << std::endl;

    utopia::USparseMatrix T =  * T_;

    utopia::disp(size(T).get(0));

    utopia::disp(size(T).get(1));

    utopia::disp(size(A_f).get(0));

    utopia::disp(size(A_m).get(0));

    utopia::set_zero_at_constraint_rows(_f_var.sys().system().get_dof_map(), A_f);

    utopia::USparseMatrix A_tot = transpose(T) * A_f * T;

    A_tot += A_m;

    utopia::disp(size(A_tot).get(0));

    utopia::disp(size(A_tot).get(1));

    utopia::UVector rhs_tot =  transpose(T) * rhs_f;

    rhs_tot += rhs_m;

    auto V_m = utopia::LibMeshFunctionSpace(utopia::make_ref(_m_var.sys().system().get_equation_systems()),_m_var.sys().system().number(), _m_var.number());

    boundary_conditions(V_m.dof_map(), A_tot, rhs_tot, rhs_m);

    utopia::disp(size(rhs_f));

    utopia::disp(size(rhs_m));

    utopia::disp(size(rhs_tot));
    
    sol_m  = utopia::local_zeros(local_size(rhs_m));

    PetscErrorCode ierr; PC _diff_problem;

    ierr = PCCreate(PETSC_COMM_WORLD, &_diff_problem); 
    CHKERRQ(ierr);
    
    ierr = PCSetType(_diff_problem,PCLU);
    CHKERRQ(ierr);
    
    ierr = PCSetOperators(_diff_problem, utopia::raw_type(A_tot),utopia::raw_type(A_tot));
    CHKERRQ(ierr);  
    
    ierr = PCFactorSetMatSolverPackage(_diff_problem,MATSOLVERMUMPS);
    CHKERRQ(ierr);
    
    ierr = PCApply(_diff_problem,utopia::raw_type(rhs_tot),utopia::raw_type(sol_m)); CHKERRQ(ierr);
    CHKERRQ(ierr);

    sol_f = T * sol_m;


    _console << "Solve_monolithic():: STOP SOLVING"  << std::endl;

    ok = true;
    
    return 1.0;

}



void
HybridModel::constraint_vec(utopia::UVector &boundary, utopia::UVector &vec)
{


    using namespace utopia;
    
    {
        Write<utopia::UVector> w_v(vec);

        Read<utopia::UVector> r_v(boundary);

            Range r = range(vec);

            for(SizeType i = r.begin(); i < r.end(); ++i) {

                if(boundary.get(i)!=0) {


                    vec.set(i, boundary.get(i));
    
                
            }
        }
    }

    synchronize(vec);
    
}

void
HybridModel::constraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat)
{


    using namespace utopia;

    typedef UTOPIA_SIZE_TYPE(UVector) SizeType;

    std::vector<SizeType> rows;
    
    {

        Read<utopia::UVector> r_v(boundary);
        
        rows.reserve(local_size(mat).get(0));

        auto r = row_range(mat);

        for(SizeType i = r.begin(); i < r.end(); ++i) {

            if(boundary.get(i)!=0) {
                rows.push_back(i);
            }
        }
    

        set_zero_rows(mat, rows, 1.);
    }    
}


void 
HybridModel::set_zero_at_constraint_rows(DofMap &dof_map, utopia::USparseMatrix &mat)
{
    bool has_constaints = true;

    using namespace utopia;

    if( dof_map.constraint_rows_begin() == dof_map.constraint_rows_end()) {
        // std::cerr << "[Warning] no zero boundary conditions to apply\n" << std::endl;
        has_constaints = false;
    }

    Size s = size(mat);
    utopia::USparseMatrix temp = mat;

    {
        Write<USparseMatrix> w_t(mat);

        each_read(temp, [&](const SizeType i, const SizeType j, const libMesh::Real value) {
            if(has_constaints && dof_map.is_constrained_dof(i)) {
                mat.set(i, j, 0.0);
            }
        });
    }
}



void
HybridModel::find_boundary(std::vector<int> &zero_rows, std::vector<int> &_dc_boundary_id){


    ConstBndNodeRange & bnd_nodes = *_fe_problem.mesh().getBoundaryNodeRange();
    NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();
    unsigned int i = 0;

    for(auto boundary = _dc_boundary_id.begin(); boundary != _dc_boundary_id.end(); ++boundary, i++)
      {
        // iterate just over boundary nodes
            for (const auto & bnode : bnd_nodes)
            {
                  libMesh::Node * current_node = bnode->_node;

                  // check if node is in active boundary list
                  if (_fe_problem.mesh().isBoundaryNode(current_node->id(), *boundary))
                  {
                    // loop over all variables at this node

                    for (auto v = 0; v < _fe_problem.getNonlinearSystemBase().nVariables(); v++)
                    {
                      const Variable & var = _nl.system().variable(v);
                      unsigned int var_num = var.number();

                      // see if this variable has any dofs at this node
                      if (current_node->n_dofs(_fe_problem.getNonlinearSystemBase().number(), var_num) > 0)
                      {
                        // check if given variable has BC on node

                        if(std::find(_dc_variables_id_m[i].begin(), _dc_variables_id_m[i].end(), var_num) != _dc_variables_id_m[i].end())
                        {

                          zero_rows.push_back(
                              current_node->dof_number(_fe_problem.getNonlinearSystemBase().number(), var_num, 0));
                        }
                    }
                }
            } 
        }
    }

    std::cout<<"Find BC with zero_rows_size: "<< zero_rows.size()<<std::endl;
}


void 
HybridModel::determine_dc_bnd_var_id(const std::vector<std::string> & BC_var){
 
     NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();

    std::vector<int> vec(_nl.nVariables());

    std::iota(vec.begin(), vec.end(), 0);

    unsigned int i;

    auto str_tmp = BC_var.begin();

    PetscFunctionBegin;
    // going over all BC_ids
    for(i = 0; str_tmp != BC_var.end(); i++, str_tmp++)
    {
        std::vector<std::string> tmp = HybridModel::split_string(*str_tmp, '-');

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
            _dc_variables_id_m.push_back(vec);
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
            _dc_variables_id_m.push_back(one_BC_id);
        }
    }

    // check if u have same number of BC_ids in both parameters
    if(_dc_variables_id_m.size() != _dc_boundary_id_m.size())
    {
        _dc_variables_id_m.clear();
        for(auto i = 0; i != _dc_boundary_id_m.size(); i++)
        {
            _dc_variables_id_m.push_back(vec);
        }
    }

    // print out what is considered for zero-ing
    std::cout<<" ------ BC CONDITIONS  ------ \n";
    unsigned int t = 0;
    //std::cout<<"_dc_variables_id.begin()"<<_dc_variables_id.size()<<std::endl;
    for(auto i = _dc_variables_id_m.begin(); i != _dc_variables_id_m.end();  t++, i++)
    {
        std::cout<<"\n BC_id:  "<< _dc_boundary_id_m[t] << "   var_ids:  ";
        std::for_each(i->begin(), i->end(), [](int i){ std::cout << i << "  " ; });
    }

}


    

     
std::vector<std::string>
HybridModel::split_string(const std::string & s, char delim)
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

 void 
 HybridModel::boundary_conditions(libMesh::DofMap &dof_map, utopia::USparseMatrix &mat, utopia::UVector &vec, utopia::UVector &rhs_values)
{ 


    determine_dc_bnd_var_id(HybridModel::split_string(_dc_var_m, ' '));
        
    find_boundary(zero_rows, _dc_boundary_id_m);

    const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

    utopia::Size ls = utopia::local_size(mat);
    utopia::Size s = utopia::size(mat);

    std::vector<int> index;

    utopia::Range rr = utopia::range(vec);


    for(int i = rr.begin(); i < rr.end(); ++i) {
        auto it = std::find(zero_rows.begin(), zero_rows.end(), i);
        if(it!=zero_rows.end()) {
            index.push_back(i);
            std::cout<<"i"<<i<<std::endl;
        }
    }

    utopia::set_zero_rows(mat, index, 1.);

    utopia::Write<utopia::UVector> w_v(vec);

    utopia::Read<utopia::UVector> r_v(rhs_values);




    utopia::Range r = range(vec);
    for(int i = r.begin(); i < r.end(); ++i) {
        auto it = std::find(zero_rows.begin(), zero_rows.end(), i);
        if(it!=zero_rows.end()) {
            auto value = rhs_values.get(i);
            vec.set(i, value);
        }
    }
    
}




void
HybridModel::CopyMatrixSolution(utopia::UVector _sol_m)
{


    FEProblemBase & _problem_m = _fe_problem;

    MooseVariable & _var_m = _problem_m.getStandardVariable(0, _m_var_name);
    
    System & _sys_m = _var_m.sys().system();
    
    MeshBase *_mesh_m = &_problem_m.mesh().getMesh();
    
    NumericVector<Number> * _solution_m = _sys_m.solution.get();
    
    PetscInt       rstart,rend;
    
    PetscScalar    tmp_sol;
    
    VecGetOwnershipRange(utopia::raw_type(_sol_m),&rstart,&rend);
    
    _sol_m = _sol_m * (-1);
    
    utopia::Read<VecT> w_d(_sol_m);

    
    {
        MeshBase::const_node_iterator it = _mesh_m->local_nodes_begin();
        const MeshBase::const_node_iterator end_it = _mesh_m->local_nodes_end();
        for ( ; it != end_it; ++it)
        {
            const Node * node = *it;
            
            for (unsigned int comp = 0;comp < node->n_comp(_sys_m.number(), _var_m.number()); comp++)
            {
                const dof_id_type from_index = node->dof_number(_sys_m.number(), _var_m.number(), comp);
                
                _solution_m->set(from_index, _sol_m.get(from_index));
                
            }
        }
    }
    
    {
        MeshBase::const_element_iterator it = _mesh_m->active_local_elements_begin();
        const MeshBase::const_element_iterator end_it = _mesh_m->active_local_elements_end();
        for ( ; it != end_it; ++it)
        {
            const Elem * elem = *it;
            for (unsigned int comp = 0; comp < elem->n_comp(_sys_m.number(), _var_m.number()); comp++)
            {
                const dof_id_type from_index = elem->dof_number(_sys_m.number(), _var_m.number(), comp);
                
                _solution_m->set(from_index, _sol_m.get(from_index));
            }
        }
    }
    
    
    _solution_m->close();
    _sys_m.update();

    Real time =  _fe_problem.dt() * _fe_problem.timeStep();

    if(_pressure)
    ExodusII_IO (*_mesh_m).write_equation_systems("matrix_p.e", _var_m.sys().system().get_equation_systems());
    else{

        if(!_ex_writer) _ex_writer = libmesh_make_unique<ExodusII_IO>(*_mesh_m);
                // A pretty update message
          libMesh::out << "\n\n*** Solving time step "
                       << _problem_m.timeStep()
                       << ", time = "
                       << time
                       << " ***"
                       << std::endl; 
        _ex_writer->write_timestep("matrix_c.e", _var_m.sys().system().get_equation_systems(),_problem_m.timeStep(), time);
    }
    
    
    
}

void
HybridModel::CopyFractureSolution(utopia::UVector _sol_f)
{
    
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);

    FEProblemBase & _f_problem =_multi_app.appProblemBase(0);

    MooseVariable & _var_f   = _f_problem.getStandardVariable(0, _f_var_name);
    
    System & _sys_f = _var_f.sys().system();

    MeshBase * _mesh_f = &_f_problem.mesh().getMesh();
    
    NumericVector<Number> * _solution_f = _sys_f.solution.get();
    
    PetscInt       rstart,rend;
    
    PetscScalar    tmp_sol;
    
    VecGetOwnershipRange(utopia::raw_type(_sol_f),&rstart,&rend);
    
    _sol_f=_sol_f*(-1);
    
    utopia::Read<VecT> w_d(_sol_f);
    
    
    {
        MeshBase::const_node_iterator it = _mesh_f->local_nodes_begin();

        const MeshBase::const_node_iterator end_it = _mesh_f->local_nodes_end();

        for ( ; it != end_it; ++it)
        {
            const Node * node = *it;
            
            for (unsigned int comp = 0;comp < node->n_comp(_sys_f.number(), _var_f.number()); comp++)
            {
                const dof_id_type to_index = node->dof_number(_sys_f.number(), _var_f.number(), comp);
                
                _solution_f->set(to_index,  _sol_f.get(to_index));
                
            }
        }
    }
    
    {
        MeshBase::const_element_iterator it = _mesh_f->active_local_elements_begin();
        const MeshBase::const_element_iterator end_it = _mesh_f->active_local_elements_end();
        for ( ; it != end_it; ++it)
        {
            const Elem * elem = *it;
            for (unsigned int comp = 0; comp < elem->n_comp(_sys_f.number(), _var_f.number()); comp++)
            {
                const dof_id_type to_index = elem->dof_number(_sys_f.number(), _var_f.number(), comp);
                _solution_f->set(to_index, _sol_f.get(to_index));
            }
        }
    }
    
    
    _solution_f->close();
    _sys_f.update();

    Real time =  _fe_problem.dt() * _fe_problem.timeStep();
    
    if(_pressure)
       ExodusII_IO (*_mesh_f).write_equation_systems("fracture_p.e", _var_f.sys().system().get_equation_systems());
    
    else{
        
       if(!_ex_writer) _ex_writer = libmesh_make_unique<ExodusII_IO>(*_mesh_f);
                // A pretty update message
          libMesh::out << "\n\n*** Solving time step "
                       << _f_problem.timeStep()
                       << ", time = "
                       << time
                       << " ***"
                       << std::endl; 
        _ex_writer->write_timestep("fracture_c.e", _var_f.sys().system().get_equation_systems(),_f_problem.timeStep(), time);
    }
    

    
}









void
HybridModel::assemble_poro_mass_matrix(FEProblemBase & _problem, std::string &_userobject_name){
    
   _console << "Assemble_Poro_Mass_matrix() begin "  << std::endl;

    utopia::USparseMatrix mass_P, lumped_mass_P;
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _problem.es().get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();    
    
    // Get a reference to our system.
    if(_fe_problem.timeStep()==1) {
        _problem.es().add_system<LinearImplicitSystem>("aux").add_variable("var",FIRST);

        _problem.es().reinit();
    }

    // Get a reference to our system.
    LinearImplicitSystem & _system = _problem.es().get_system<LinearImplicitSystem>("aux");
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);
    
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    
    SparseMatrix<Number> & matrix_M_p = *_system.matrix;
    
    //_console << "is matrix closed: " << matrix_A.closed() << std::endl;
    
    // The element mass matrix.
    DenseMatrix<Number> Me_p;
    
    // A  Gauss quadrature rule for numerical integration.
    QGauss qrule (dim, fe_type.default_quadrature_order());
    
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);
    
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();
    
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    
    const DofMap & dof_map = _system.get_dof_map();
    
    std::vector<dof_id_type> dof_indices;
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    // first we need to manually zero the matrix
    matrix_M_p.zero();
    
    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;

        Elem * ele = *el;
        
        fe->reinit (elem);

        dof_map.dof_indices(elem, dof_indices);
        
        Me_p.resize (dof_indices.size(), dof_indices.size());
        
        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                        Me_p(i,j) +=   ComputeMaterialProprties(elem) * JxW[qp] * phi[i][qp] * phi[j][qp];
            }
        }
    }

        
        dof_map.constrain_element_matrix(Me_p,dof_indices,false);

        matrix_M_p.add_matrix(Me_p, dof_indices);

        
    }


    matrix_M_p.close();

    _console << "Assemble_Poro_Mass_matrix() begin 2"  << std::endl;

    utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();
        
    utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

    utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M_p), mass_temp);

    mass_P = mass_temp;

    lumped_mass_P = diag(sum(mass_P,1));

    _console << "Assemble_Poro_Mass_matrix() end "  << std::endl;

    const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name)).setOperator() = std::make_shared<utopia::USparseMatrix>(lumped_mass_P);

    
}



void
HybridModel::stabilize_matrix(utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix)
{ 
    

   _console << "Stabilize A matrix:: begin  "  << std::endl;

    utopia::USparseMatrix A_0_t;
    
    A_0_t = utopia::transpose(A_0);

    S_matrix = A_0_t;

    S_matrix*=0;
    
    
    {
        utopia::Read<utopia::USparseMatrix>  r_s(A_0), r_s_t(A_0_t);
        utopia::Write<utopia::USparseMatrix> w_s(S_matrix);
        utopia::each_read(A_0, [&](const utopia::SizeType i, const utopia::SizeType j, double value){
            if(i!=j)
            {
                double value_1 = 1.0 * value;      
                
                double value_2 = 1.0 * A_0_t.get(i,j);

                Real max=std::max(value_1,value_2);

                if (max>0.0){

                    max*=-1.0;

                    S_matrix.set(i,j,max);
                }
            }
                
           
        });
    }
    
    
    
    
    utopia::UVector diag_elem = -1.0 * sum(S_matrix,1);
    
    utopia::USparseMatrix S_diag=diag(diag_elem);

    
    S_matrix+=S_diag;
    
    _console << "Stabilize A matrix:: end  "  << std::endl;




  }
Real
HybridModel::ComputeMaterialProprties(const Elem *elem){
    
   // _console << "_vector_p.size()"  << _vector_p.size() <<std::endl;

Real poro=0.0;

    for(int ll=0; ll<_vector_p.size(); ll++){
        if (elem->subdomain_id()==_vector_p[ll]) {

            poro = _vector_value[ll]; 
            //std::cout<<"poro"<<poro<<std::endl;
        }
    }
  
    return poro;   
}





int 
HybridModel::solve_transport_monolithic(bool & ok)
{        
        

    using namespace utopia;
    
   

    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
        
    FEProblemBase & _m_problem = _multi_app.problemBase();
    
    FEProblemBase & _f_problem = _multi_app.appProblemBase(0);

    auto &_m_sys = _m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    NonlinearSystemBase & _nl_m = _m_problem.getNonlinearSystemBase();
    
    auto &_f_sys = _f_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    NonlinearSystemBase & _nl_f = _f_problem.getNonlinearSystemBase();



    if  (_fe_problem.timeStep()==1)
     {
              
        _console << "solve_transport_stabilize_static_condens::time_step "  << _fe_problem.timeStep() <<std::endl;     
        
        
        libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_m_sys.matrix);

        libMesh::PetscMatrix<libMesh::Number> *petsc_mat_f = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_f_sys.matrix);

        _m_problem.computeJacobianSys(_m_sys, *_nl_m.currentSolution(), *petsc_mat_m);

        USparseMatrix A_m_t, A_f_t;

        UVector rhs_f_t, rhs_m_t, rhs_m_c, rhs_f_c, c_m_old, c_f_old;


        utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_m).mat(), A_m_t);

        _m_problem.computeResidualSys(_m_sys, *_nl_m.currentSolution(), _nl_m.RHS());

        utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_m.RHS()), rhs_m_t);   

        _f_problem.computeJacobianSys(_f_sys, *_nl_f.currentSolution(), *petsc_mat_f);

        utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_f).mat(), A_f_t);

        _f_problem.computeResidualSys(_f_sys, *_nl_f.currentSolution(), _nl_f.RHS());

        utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_f.RHS()), rhs_f_t);



        rhs_m_c = rhs_m_t; 

        rhs_f_c = rhs_f_t; 



        Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();


        double inv_dt = 1.0/dt;

        c_m_old  = utopia::local_zeros(local_size(rhs_m_t));
        c_f_old  = utopia::local_zeros(local_size(rhs_f_t));

        UVector  mass_c_m_old, mass_c_f_old, mass_c_m_old_dot, mass_c_f_old_dot;

        USparseMatrix mass_lumped_mp = *const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_M_m)).getOperator();
        USparseMatrix mass_lumped_mpc =  mass_lumped_mp;        
        constraint_mat(rhs_m_c,  mass_lumped_mpc);       
        mass_c_m_old = mass_lumped_mp * c_m_old;
        mass_c_m_old_dot = mass_c_m_old * inv_dt;
        rhs_m_t = - 1.0 * mass_c_m_old_dot;

        
        USparseMatrix mass_lumped_fp = *const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_M_f)).getOperator();
        USparseMatrix mass_lumped_fpc =  mass_lumped_fp;     
        constraint_mat(rhs_f_c,  mass_lumped_fpc);
        mass_c_f_old = mass_lumped_fp * c_f_old;
        mass_c_f_old_dot = mass_c_f_old * inv_dt;
        rhs_f_t = - 1.0 * mass_c_f_old_dot;

        _console << "Solve_monolithic():: START SOLVING"  << std::endl;

        utopia::USparseMatrix T =  * T_;

        utopia::USparseMatrix A_tot = transpose(T) * A_f_t * T;

        A_tot += A_m_t;

        utopia::USparseMatrix   S_m_cons = A_m_t;

        S_m_cons*=0;

        stabilize_matrix(A_tot, S_m_cons);

        USparseMatrix A_sc = mass_lumped_mp * inv_dt + transpose(T) * mass_lumped_fp * inv_dt * T + A_tot + S_m_cons;
        
        UVector rhs_sc     = rhs_m_t + transpose(T) * rhs_f_t;

        constraint_mat(rhs_m_c,  A_sc);

        constraint_vec(rhs_m_c,  rhs_sc);

        _A_t_store   = std::make_shared<USparseMatrix>(A_sc);

        const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_A_tot)).setOperator();


        sol_m  = utopia::local_zeros(local_size(rhs_m_t));


        _ksp_ptr = (KSP_PARROT *)ksp->data;

        int factorized=_ksp_ptr[0].factorized[0];

        std::cout<<factorized<<std::endl;

        PetscErrorCode ierr; 
       
        
        ierr = KSPGetOperators(ksp,&raw_type(A_sc),&raw_type(A_sc));CHKERRQ(ierr);
        
        if (factorized==0)
        {
            _ksp_ptr[0].factorized[0]=1;
    
        
            PetscErrorCode ierr;
            
            PCSetType(_ksp_ptr[0].local_pc[0],PCLU);
            PCSetOperators(_ksp_ptr[0].local_pc[0],utopia::raw_type(A_sc),utopia::raw_type(A_sc));
            PCFactorSetMatSolverPackage(_ksp_ptr[0].local_pc[0],MATSOLVERMUMPS); //MATSOLVERSUPERLU_DIST);MATSOLVERMUMPS
                
            std::cout<<"start factorizing?\n";
            auto t_start = std::chrono::high_resolution_clock::now();
            PCSetUp(_ksp_ptr[0].local_pc[0]);
            auto t_end = std::chrono::high_resolution_clock::now();
            std::cout<<"done factorizing?\n";
                std::cout<<"fact time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";
        }
        
        PCSetReusePreconditioner(_ksp_ptr[0].local_pc[0], PETSC_TRUE);
        std::cout<<"start solving?\n";
        auto t_start = std::chrono::high_resolution_clock::now();
        PCApply(_ksp_ptr[0].local_pc[0],utopia::raw_type(rhs_sc),utopia::raw_type(sol_m));
        auto t_end = std::chrono::high_resolution_clock::now();
        std::cout<<"done solving?\n";
        std::cout<<"solve time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";

        _ksp_ptr[0].local_pc=NULL;

        Vec r;
        VecDuplicate(utopia::raw_type(rhs_sc),&r);
        MatResidual(utopia::raw_type(A_sc),utopia::raw_type(rhs_sc),utopia::raw_type(sol_m),r);
        PetscReal norm;
        VecNorm(r,NORM_2,&norm);
        std::cout<<"qui "<<norm<<std::endl;
        PetscPrintf(PETSC_COMM_WORLD,"   %14.12e \n", norm);

        ksp->its    = 1;
        ksp->reason = KSP_CONVERGED_ITS;



        //auto op = std::make_shared<utopia::Factorization<utopia::USparseMatrix, utopia::UVector> >(MATSOLVERMUMPS,PCLU);
        //const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_A_tot)).getVoidPointer() = op;
        // PetscErrorCode ierr; PC _diff_problem;

        // ierr = PCCreate(PETSC_COMM_WORLD, &_diff_problem); 
        // CHKERRQ(ierr);
        
        // ierr = PCSetType(_diff_problem,PCLU);
        // CHKERRQ(ierr);
        
        // ierr = PCSetOperators(_diff_problem, utopia::raw_type(A_tot),utopia::raw_type(A_sc));
        // CHKERRQ(ierr);  
        
        // ierr = PCFactorSetMatSolverPackage(_diff_problem,MATSOLVERMUMPS);
        // CHKERRQ(ierr);
        
        // ierr = PCApply(_diff_problem,utopia::raw_type(rhs_sc),utopia::raw_type(sol_m)); CHKERRQ(ierr);
        // CHKERRQ(ierr);

        // op->update(_A_t_store);

        // op->apply(rhs_sc, sol_m);

        sol_f = T * sol_m;

        sol_m *= 1;

        sol_f = T * sol_m;  


    
    }




  
        _console << "solve_transport_stabilize_static_condens::time_step "  << _fe_problem.timeStep() <<std::endl;     
        
        auto sol_m_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_m_sys.old_local_solution.get())->vec();

        auto sol_f_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_f_sys.old_local_solution.get())->vec();

        UVector rhs_f_t, rhs_m_t, rhs_m_c, rhs_f_c, c_m_old, c_f_old;
        UVector  mass_c_m_old, mass_c_f_old, mass_c_m_old_dot, mass_c_f_old_dot;

        utopia::convert(sol_m_old, c_m_old);

        utopia::convert(sol_f_old, c_f_old);


        Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();
        double inv_dt = 1.0/dt;

        

        USparseMatrix mass_lumped_mp = *const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_M_m)).getOperator();
        USparseMatrix mass_lumped_mpc =  mass_lumped_mp;        
        constraint_mat(rhs_m_c,  mass_lumped_mpc);       
        mass_c_m_old = mass_lumped_mp * c_m_old;
        mass_c_m_old_dot = mass_c_m_old * inv_dt;
        rhs_m_t = - 1.0 * mass_c_m_old_dot;

        
        USparseMatrix mass_lumped_fp = *const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_M_f)).getOperator();
        USparseMatrix mass_lumped_fpc =  mass_lumped_fp;     
        constraint_mat(rhs_f_c,  mass_lumped_fpc);
        mass_c_f_old = mass_lumped_fp * c_f_old;
        mass_c_f_old_dot = mass_c_f_old * inv_dt;
        rhs_f_t = - 1.0 * mass_c_f_old_dot;


        _A_t_store = const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_A_tot)).getOperator();

        utopia::USparseMatrix T =  * T_;

        UVector rhs_sc     = rhs_m_t + transpose(T) * rhs_f_t;

        utopia::USparseMatrix A_sc = * _A_t_store;


        PCSetReusePreconditioner(_ksp_ptr[0].local_pc[0], PETSC_TRUE);
        std::cout<<"start solving?\n";
        auto t_start = std::chrono::high_resolution_clock::now();
        PCApply(_ksp_ptr[0].local_pc[0],utopia::raw_type(rhs_sc),utopia::raw_type(sol_m));
        auto t_end = std::chrono::high_resolution_clock::now();
        std::cout<<"done solving?\n";
        std::cout<<"solve time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";

        _ksp_ptr[0].local_pc=NULL;

        Vec r;
        VecDuplicate(utopia::raw_type(rhs_sc),&r);
        MatResidual(utopia::raw_type(A_sc),utopia::raw_type(rhs_sc),utopia::raw_type(sol_m),r);
        PetscReal norm;
        VecNorm(r,NORM_2,&norm);
        std::cout<<"qui "<<norm<<std::endl;
        PetscPrintf(PETSC_COMM_WORLD,"   %14.12e \n", norm);

        ksp->its    = 1;
        ksp->reason = KSP_CONVERGED_ITS;

        //sol_m  = utopia::local_zeros(local_size(rhs_m_t));

        //auto op = std::static_pointer_cast< Factorization<USparseMatrix, UVector> >(const_cast<StoreUtopiaOperators&>(_fe_problem.getUserObjectTempl<StoreUtopiaOperators>(_userobject_name_A_tot)).getVoidPointer());

        // PetscErrorCode ierr; PC _diff_problem;

        // ierr = PCCreate(PETSC_COMM_WORLD, &_diff_problem); 
        // CHKERRQ(ierr);
        
        // ierr = PCSetType(_diff_problem,PCLU);
        // CHKERRQ(ierr);
        
        // ierr = PCSetOperators(_diff_problem, utopia::raw_type(A_sc),utopia::raw_type(A_sc));
        // CHKERRQ(ierr);  
        
        // ierr = PCFactorSetMatSolverPackage(_diff_problem,MATSOLVERMUMPS);
        // CHKERRQ(ierr);
        
        // ierr = PCApply(_diff_problem,utopia::raw_type(rhs_sc),utopia::raw_type(sol_m)); CHKERRQ(ierr);
        // CHKERRQ(ierr);
        
        //op->apply(rhs_sc, sol_m);

        sol_f = T * sol_m;

        sol_m *= 1;

        sol_f = T * sol_m;        
        
        ok = true;

        return 1.0;
    
}











  



