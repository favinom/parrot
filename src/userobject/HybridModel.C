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




using namespace std;

registerMooseObject("parrotApp",HybridModel);


typedef utopia::USparseMatrix SparseMatT;
typedef utopia::UVector VecT;




template <>
InputParameters
validParams<HybridModel>()
{
    InputParameters params = validParams<GeneralUserObject>();
    
    // // params.addRequiredParam<std::vector<int>>("block_id",
    // //                                                  "The name of the nodeset to create");
    // // params.addRequiredParam<std::vector<Real>>("value_p",
    // //                                                  "The name of the nodeset to create");
    // // params.addRequiredParam<AuxVariableName>("lagrange_variable",
    // //                                          "The auxiliary variable to store the transferred values in.");
    params.addRequiredParam<VariableName>("matrix_variable",
                                          "The variable to transfer from.");
    params.addRequiredParam<VariableName>("fracture_variable",
                                          "The variable to transfer from.");
    // params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    // // params.addParam<bool>("pressure","false","put true if you solve pressure system");
    // // params.addParam<bool>("transport","false","put true if you transport system");
    // params.addParam<bool>("biorth","false","put true for biorth bases");
    // // params.addParam<bool>("stabilize","false","put true for FCT stabilizing");
    // // params.addParam<bool>("constraint_m","false","put true if matrix has Dirichlet BC");
    // // params.addParam<bool>("constraint_f","false","put true if fibres has Dirichlet BC");
    // params.addRequiredParam<MultiAppName>("multi_app", "The MultiApp's name in your input file!");
    
    return params;
}

HybridModel::HybridModel(const InputParameters & parameters) :
GeneralUserObject(parameters),
 _f_var_name(getParam<VariableName>("fracture_variable")),
 _m_var_name(getParam<VariableName>("matrix_variable")),
_operator_storage(getUserObject<StoreOperators>("operator_userobject")),
// //_pressure(getParam<bool>("pressure")),
// //_transport(getParam<bool>("transport")),
// _biorth(getParam<bool>("biorth")),
// // _vector_p(getParam<std::vector<int>>("block_id")),
// // _vector_value(getParam<std::vector<Real>>("value_p")),
_multiapp_name(getParam<MultiAppName>("multi_app"))
{
    
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

    // if(!m_problem.hasUserObject(_userobject_name_S)){

    //         std::string class_name = "StoreTransferOperators";
            
    //         auto params = m_problem.getMooseApp().getFactory().getValidParams(class_name);
            
    //         params.set<bool>("use_displaced_mesh") = false;
            
    //         params.set<ExecFlagEnum>("execute_on") = "initial";

    //         // m_problem.addUserObject("StoreOperators", _userobject_name_S, params);

    //         // m_problem.addUserObject("StoreOperators", _userobject_name_A_m, params);

    //         // m_problem.addUserObject("StoreOperators", _userobject_name_A_f, params);

    //         // m_problem.addUserObject("StoreOperators", _userobject_name_T, params);

    // }

    

        if  (m_problem.timeStep()==1){

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

                utopia::set_zero_at_constraint_rows(_m_var.sys().system().get_dof_map(), *T_);
        }


    // if (_transport){

    //     assemble_mass_matrix(_porosity_m, m_problem,  mass_m, mass_lumped_m);

    //     assemble_mass_matrix(_porosity_f, f_problem, mass_f, mass_lumped_f);

    //     assemble_poro_mass_matrix(_porosity_m, m_problem,  mass_mp, mass_lumped_mp);

    //     assemble_poro_mass_matrix(_porosity_f, f_problem, mass_f, mass_lumped_fp);
    // }


}



bool
HybridModel::solve(){


    solve_pressure_monolithic();

    return ok;
    
}


void
HybridModel::execute()
{

    solve(); 
        
    if (ok){
        
        // CopyMatrixSolution(sol_m);
        
        // CopyFractureSolution(sol_f);
    }


   
}


bool HybridModel::solve_pressure_monolithic()
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

    utopia::USparseMatrix A_tot = transpose(T) * A_f * T + A_m;

    utopia::UVector rhs_tot     = rhs_m + transpose(T) * rhs_f;
    
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
    return true;

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

    //if(_pressure)
       ExodusII_IO (*_mesh_m).write_equation_systems("matrix_p.e", _var_m.sys().system().get_equation_systems());
    // else
    //    ExodusII_IO (*_mesh_m).write_equation_systems("matrix_c.e", _var_m.sys().system().get_equation_systems());
    
    
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
    
    //if(_pressure)
       ExodusII_IO (*_mesh_f).write_equation_systems("fracture_p.e", _var_f.sys().system().get_equation_systems());
    
    // else{
        
    //    if(!_ex_writer) _ex_writer = libmesh_make_unique<ExodusII_IO>(*_mesh_f);
    //             // A pretty update message
    //       libMesh::out << "\n\n*** Solving time step "
    //                    << _f_problem.timeStep()
    //                    << ", time = "
    //                    << time
    //                    << " ***"
    //                    << std::endl; 
    //     _ex_writer->write_timestep("fracture_c.e", _var_f.sys().system().get_equation_systems(),_f_problem.timeStep(), time);
    //}
    

    
}



// void
// FractureAppNConforming::constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
// {


//     using namespace utopia;
    
//     {
//         Write<utopia::UVector> w_v(vec);

//         Read<utopia::UVector> r_v(boundary);

//         if(has_constaints) {

//             Range r = range(vec);

//             for(SizeType i = r.begin(); i < r.end(); ++i) {

   

//                 if(boundary.get(i)!=0) {


//                     vec.set(i, boundary.get(i));
    
//                 }
//             }
//         }
//     }

//     synchronize(vec);
    
// }


// void
// FractureAppNConforming::unconstraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
// {


//     using namespace utopia;
    
//     {
//         Write<utopia::UVector> w_v(vec);

//         Read<utopia::UVector> r_v(boundary);

//         if(has_constaints) {

//             Range r = range(vec);

//             for(SizeType i = r.begin(); i < r.end(); ++i) {

   

//                 if(boundary.get(i)!=0) {


//                     vec.set(i, 0.0);
    
//                 }
//             }
//         }
//     }

//     synchronize(vec);
    
// }

// void
// FractureAppNConforming::one_constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
// {


//     using namespace utopia;
    
//     {
//         Write<utopia::UVector> w_v(vec);

//         Read<utopia::UVector> r_v(boundary);

//         if(has_constaints) {

//             Range r = range(vec);

//             for(SizeType i = r.begin(); i < r.end(); ++i) {

   

//                 if(boundary.get(i)!=0) {


//                     vec.set(i, 1.0);
    
//                 }
//             }
//         }
//     }

//     synchronize(vec);
    
// }


// void
// FractureAppNConforming::constraint_concentration_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints)
// {
// //    auto &V = space->space().last_subspace();

//     using namespace utopia;

//     typedef UTOPIA_SIZE_TYPE(UVector) SizeType;

//     std::vector<SizeType> rows;
    
//     {

//         Read<utopia::UVector> r_v(boundary);
        
//         rows.reserve(local_size(mat).get(0));

//         if(has_constaints) {

//             auto r = row_range(mat);

//             for(SizeType i = r.begin(); i < r.end(); ++i) {

//                 //std::cout<<"value "<< boundary.get(i) <<std::endl;

//                 if(boundary.get(i)!=0) {

//                     //std::cout<<"i "<< i <<"valpos->second "<<boundary.get(i)<<std::endl;

//                     // if(valpos != rhs_values.end()) {
//                     rows.push_back(i);
//                     // }
//                 }
//             }
//         }

//         set_zero_rows(mat, rows, 1.);
//     }    
// }




// void
// FractureAppNConforming::assemble_mass_matrix(double porosity, FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix, utopia::USparseMatrix &lumped_mass_matrix){
    
//    _console << "Assemble_Mass_matrix() begin "  << std::endl;
    
//     // Get a constant reference to the mesh object.
//     const MeshBase & mesh = _problem.es().get_mesh();
    
//     // The dimension that we are running.
//     const unsigned int dim = mesh.mesh_dimension();    
    
//     // Get a reference to our system.
//     auto & _system = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
//     // Get a constant reference to the Finite Element type
//     // for the first (and only) variable in the system.
//     FEType fe_type = _system.get_dof_map().variable_type(0);
    
//     UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    
//     SparseMatrix<Number> & matrix_M = *_system.matrix;
    
//     // The element mass matrix.
//     DenseMatrix<Number> Me;

//     //DenseMatrix<Number> Me_p;
    
//     // A  Gauss quadrature rule for numerical integration.
//     QGauss qrule (dim, fe_type.default_quadrature_order());
    
//     // Tell the finite element object to use our quadrature rule.
//     fe->attach_quadrature_rule (&qrule);
    
//     // The element Jacobian * quadrature weight at each integration point.
//     const std::vector<Real> & JxW = fe->get_JxW();
    
//     // The element shape functions evaluated at the quadrature points.
//     const std::vector<std::vector<Real> > & phi = fe->get_phi();
    
//     const DofMap & dof_map = _problem.getNonlinearSystemBase().dofMap();
    
//     std::vector<dof_id_type> dof_indices;
    
//     MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

//     const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
//     matrix_M.zero();
    
//     for ( ; el != end_el; ++el)
//     {
//         const Elem * elem = *el;

//         Elem * ele = *el;
        
//         fe->reinit (elem);

//         dof_map.dof_indices(elem, dof_indices);
        
//         Me.resize (dof_indices.size(), dof_indices.size());

//         // Me_p.resize (dof_indices.size(), dof_indices.size());
        
//         for (unsigned int qp=0; qp<qrule.n_points(); qp++){

//             for (unsigned int i=0; i<phi.size(); i++){

//                 for (unsigned int j=0; j<phi.size(); j++){

//                         Me(i,j) += JxW[qp] * phi[i][qp] * phi[j][qp];
//                 }      
//             }
//         }

        
//         dof_map.constrain_element_matrix(Me,dof_indices,false);

//         matrix_M.add_matrix (Me, dof_indices);
        
//     }


    
//     matrix_M.close();

//     utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();
        
//     utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

//     utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M), mass_temp);

//     mass_matrix =  mass_temp;

//     lumped_mass_matrix = diag(sum(mass_matrix,1));


//     _console << "Assemble_Mass_matrix() end "  << std::endl;

    
// }


// void
// FractureAppNConforming::assemble_poro_mass_matrix(double porosity, FEProblemBase & _problem, utopia::USparseMatrix &mass_matrixp, utopia::USparseMatrix &lumped_mass_matrixp){
    
//    _console << "Assemble_Poro_Mass_matrix() begin "  << std::endl;
    
//     // Get a constant reference to the mesh object.
//     const MeshBase & mesh = _problem.es().get_mesh();
    
//     // The dimension that we are running.
//     const unsigned int dim = mesh.mesh_dimension();    
    
//     // Get a reference to our system.
//     if(_fe_problem.timeStep()==1) {
//         _problem.es().add_system<LinearImplicitSystem>("aux").add_variable("var",FIRST);

//         _problem.es().reinit();
//     }

//     // Get a reference to our system.
//     LinearImplicitSystem & _system = _problem.es().get_system<LinearImplicitSystem>("aux");
    
//     // Get a constant reference to the Finite Element type
//     // for the first (and only) variable in the system.
//     FEType fe_type = _system.get_dof_map().variable_type(0);
    
//     UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    
//     SparseMatrix<Number> & matrix_M_p = *_system.matrix;
    
//     //_console << "is matrix closed: " << matrix_A.closed() << std::endl;
    
//     // The element mass matrix.
//     DenseMatrix<Number> Me_p;
    
//     // A  Gauss quadrature rule for numerical integration.
//     QGauss qrule (dim, fe_type.default_quadrature_order());
    
//     // Tell the finite element object to use our quadrature rule.
//     fe->attach_quadrature_rule (&qrule);
    
//     // The element Jacobian * quadrature weight at each integration point.
//     const std::vector<Real> & JxW = fe->get_JxW();
    
//     // The element shape functions evaluated at the quadrature points.
//     const std::vector<std::vector<Real> > & phi = fe->get_phi();
    
//     const DofMap & dof_map = _system.get_dof_map();
    
//     std::vector<dof_id_type> dof_indices;
    
//     MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

//     const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
//     // first we need to manually zero the matrix
//     matrix_M_p.zero();
    
//     for ( ; el != end_el; ++el)
//     {
//         const Elem * elem = *el;

//         Elem * ele = *el;
        
//         fe->reinit (elem);

//         dof_map.dof_indices(elem, dof_indices);
        
//         Me_p.resize (dof_indices.size(), dof_indices.size());
        
//         for (unsigned int qp=0; qp<qrule.n_points(); qp++){

//             for (unsigned int i=0; i<phi.size(); i++){

//                 for (unsigned int j=0; j<phi.size(); j++){

//                         Me_p(i,j) +=   ComputeMaterialProprties(elem) * JxW[qp] * phi[i][qp] * phi[j][qp];
//             }
//         }
//     }

        
//         dof_map.constrain_element_matrix(Me_p,dof_indices,false);

//         matrix_M_p.add_matrix(Me_p, dof_indices);

        
//     }


//     matrix_M_p.close();

//     _console << "Assemble_Poro_Mass_matrix() begin 2"  << std::endl;

//     utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();
        
//     utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

//     utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M_p), mass_temp);

//     mass_matrixp = mass_temp;

//     lumped_mass_matrixp = diag(sum(mass_matrixp,1));

//     _console << "Assemble_Poro_Mass_matrix() end "  << std::endl;

    
// }



// void
// FractureAppNConforming::stabilize_matrix(utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix)
// { 
    

//    _console << "Stabilize A matrix:: begin  "  << std::endl;

//     utopia::USparseMatrix A_0_t;
    
//     A_0_t = utopia::transpose(A_0);

//     S_matrix = A_0_t;

//     S_matrix*=0;
    
    
//     {
//         utopia::Read<utopia::USparseMatrix>  r_s(A_0), r_s_t(A_0_t);
//         utopia::Write<utopia::USparseMatrix> w_s(S_matrix);
//         utopia::each_read(A_0, [&](const utopia::SizeType i, const utopia::SizeType j, double value){
//             if(i!=j)
//             {
//                 double value_1 = 1.0 * value;      
                
//                 double value_2 = 1.0 * A_0_t.get(i,j);

//                 Real max=std::max(value_1,value_2);

//                 if (max>0.0){

//                     max*=-1.0;

//                     S_matrix.set(i,j,max);
//                 }
//             }
                
           
//         });
//     }
    
    
    
    
//     utopia::UVector diag_elem = -1.0 * sum(S_matrix,1);
    
//     utopia::USparseMatrix S_diag=diag(diag_elem);

    
//     S_matrix+=S_diag;
    
//     _console << "Stabilize A matrix:: end  "  << std::endl;




//   }





//  void FractureAppNConforming::compute_static_condenstation_transport()
// {        
        

//     using namespace utopia;
    
   

//     MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
        
//     FEProblemBase & _m_problem = _multi_app.problemBase();
    
//     FEProblemBase & _f_problem = _multi_app.appProblemBase(0);

//     auto &_m_sys = _m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

//     NonlinearSystemBase & _nl_m = _m_problem.getNonlinearSystemBase();
    
//     auto &_f_sys = _f_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

//     NonlinearSystemBase & _nl_f = _f_problem.getNonlinearSystemBase();



//     if  (_fe_problem.timeStep()==1)
//      {
              
//         _console << "solve_transport_stabilize_static_condens::time_step "  << _fe_problem.timeStep() <<std::endl;     
        
//         // Pointer to underlying PetscMatrix type
//         libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_m_sys.matrix);

//         libMesh::PetscMatrix<libMesh::Number> *petsc_mat_f = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_f_sys.matrix);



//         _m_problem.computeJacobianSys(_m_sys, *_nl_m.currentSolution(), *petsc_mat_m);

//         utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_m).mat(), A_m_t);

//         _m_problem.computeResidualSys(_m_sys, *_nl_m.currentSolution(), _nl_m.RHS());

//         utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_m.RHS()), rhs_m_t);

        

//         _f_problem.computeJacobianSys(_f_sys, *_nl_f.currentSolution(), *petsc_mat_f);

//         utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_f).mat(), A_f_t);

//         _f_problem.computeResidualSys(_f_sys, *_nl_f.currentSolution(), _nl_f.RHS());

//         utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_f.RHS()), rhs_f_t);



//         rhs_m_c = rhs_m_t; 

//         rhs_f_c = rhs_f_t; 



//         Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();


//         double inv_dt = 1.0/dt;


//         USparseMatrix S_m = A_m_t;

//         S_m*=0;


//         USparseMatrix S_f = A_f_t;

//         S_f*=0;


//         c_m  = utopia::local_zeros(local_size(rhs_m_t));

//         c_f  = utopia::local_zeros(local_size(rhs_f_t));

//         UVector  mass_c_m_old, mass_c_f_old;

//         USparseMatrix mass_lumped_mpc =  mass_lumped_mp;
        
//         constraint_concentration_mat(rhs_m_c,  mass_lumped_mpc, _constraint_m);
        

//         mass_c_m_old = mass_lumped_mp * c_m;

//         mass_c_m_old_dot = mass_c_m_old * inv_dt;

//         rhs_m_t = - 1.0 * mass_c_m_old_dot;
    



//         USparseMatrix mass_lumped_fpc =  mass_lumped_fp;
        
//         constraint_concentration_mat(rhs_f_c,  mass_lumped_fpc, _constraint_f);


//         mass_c_f_old = mass_lumped_fp * c_f;

//         mass_c_f_old_dot = mass_c_f_old * inv_dt;

//         rhs_f_t = - 1.0 * mass_c_f_old_dot;
//         // constraint_concentration_vec(rhs_f_c,  rhs_f_t, _constraint_f);



//         const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_2)).setTransferOperator()= std::make_shared<USparseMatrix>(S_f);


       
//         auto op = std::make_shared<utopia::Factorization<utopia::USparseMatrix, utopia::UVector> >(MATSOLVERMUMPS,PCLU);

//         USparseMatrix D_b = diag(sum(B, 1));

//         UVector diag_elem_b = 1./sum((B),1);
        
//         USparseMatrix Dinv = diag(diag_elem_b);

//         USparseMatrix T =  Dinv * B;

//         const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_T)).setTransferOperator()= std::make_shared<USparseMatrix>(T);

//         USparseMatrix A_m_tot_cons = transpose(T) * A_f_t * T + A_m_t;

//         USparseMatrix S_m_cons = A_m_t;

//         S_m_cons*=0;

//         stabilize_matrix(A_m_tot_cons, S_m_cons);


//         const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).setTransferOperator()= std::make_shared<USparseMatrix>(S_m_cons);


//         USparseMatrix A_sc = mass_lumped_mp * inv_dt + transpose(T) * mass_lumped_fp * inv_dt * T + A_m_tot_cons + S_m_cons;
        
//         UVector rhs_sc     = rhs_m_t + transpose(T) * rhs_f_t;

//         constraint_concentration_mat(rhs_m_c,  A_sc, _constraint_m);

//         constraint_concentration_vec(rhs_m_c,  rhs_sc, _constraint_m);

//         _A_t_store   = std::make_shared<USparseMatrix>(A_sc);

//         const_cast<StoreTransferOperators&>(_operator_storage).setTransferOperator() = _A_t_store;

//         const_cast<StoreTransferOperators&>(_operator_storage).getVoidPointer() = op;

//         dlagr_t = utopia::local_zeros(local_size(B).get(0));   

//         lagr_t = utopia::local_zeros(local_size(dlagr_t));

//         c_m = utopia::local_zeros(local_size(rhs_m_t));

//         c_f = utopia::local_zeros(local_size(rhs_f_t));

//         op->update(_A_t_store);
        
//         op->apply(rhs_sc, c_m);

//         c_m*=1;

//         c_f = T * c_m;

//         CopyMatrixSolution(c_m);
        
//         CopyFractureSolution(c_f);

        
    
//     }


//     if (_fe_problem.timeStep()>1)
//     {


//          _console << "solve_transport_stabilize_static_condens::time_step "  << _fe_problem.timeStep() <<std::endl;     
        
//         auto sol_m_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_m_sys.old_local_solution.get())->vec();

//         auto sol_f_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_f_sys.old_local_solution.get())->vec();

//         utopia::convert(sol_m_old, c_m_old);

//         utopia::convert(sol_f_old, c_f_old);



//         // rhs_m_c = rhs_m_t; 

//         // rhs_f_c = rhs_f_t; 

//         // disp("c_m_old");
        
//         // disp(c_m_old);

//         // disp("rhs_m_t");
        
//         // disp(rhs_m_t);

//         Real dt = static_cast<Transient*>(_m_problem.getMooseApp().getExecutioner())->getDT();

//         _console << "transport_monolithic again"  << std::endl;

//         double inv_dt = 1.0/dt;

//         UVector  mass_c_m_old,  mass_c_f_old;

//         mass_c_m_old = mass_lumped_mp * c_m_old;

//         mass_c_f_old = mass_lumped_fp * c_f_old;

//         mass_c_f_old_dot = mass_c_f_old * inv_dt;

//         mass_c_m_old_dot = mass_c_m_old * inv_dt;

//         rhs_m_t = -1.0 * mass_c_m_old_dot;

//         rhs_f_t = -1.0 * mass_c_f_old_dot;

   

//         //constraint_concentration_vec(rhs_m_c,  rhs_m_t, _constraint_m);
    
//         // constraint_concentration_vec(rhs_f_c,  rhs_f_t, _constraint_f);



        
//         auto op = std::static_pointer_cast< Factorization<USparseMatrix, UVector> >(const_cast<StoreTransferOperators&>(_operator_storage).getVoidPointer());

        
//         dlagr_t = utopia::local_zeros(local_size(B).get(0));   

//         lagr_t = utopia::local_zeros(local_size(dlagr_t));

//         utopia::UVector rhs_t = utopia::blocks(rhs_m_t, rhs_f_t, dlagr_t);
        
//         utopia::UVector sol_t = blocks(c_m, c_f, lagr_t);

//         _A_t_store = const_cast<StoreTransferOperators&>(_operator_storage).getTransferOperator();

//         c_m = utopia::local_zeros(local_size(rhs_m_t));

//         c_f = utopia::local_zeros(local_size(rhs_m_t));

//         USparseMatrix T = *const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_T)).getTransferOperator();

//         UVector rhs_sc  = rhs_m_t + transpose(T) * rhs_f_t;

//         //op->atol(1e-4);
//         //op->rtol(1e-3);

//         op->update(_A_t_store);

//         constraint_concentration_vec(rhs_m_c,  rhs_sc, _constraint_m);

//         c_m = utopia::local_zeros(local_size(rhs_m_t));

//         c_f = utopia::local_zeros(local_size(rhs_f_t));
        
//         op->apply(rhs_sc, c_m);

     



//         //utopia::undo_blocks(sol_t, c_m, c_f, lagr_t);

//         c_m*=1;

//         c_f = T * c_m;



//         // CopyMatrixSolution(c_m);
        
//         // CopyFractureSolution(c_f);
        
//         // op->update(_A_t_store);
        
//         // op->apply(rhs_t, sol_t);

//         // utopia::undo_blocks(sol_t, c_m, c_f, lagr_t);

//         // c_m*=1;

//         // c_f*=1;


//         utopia::UVector rhs_m_dot, rhs_f_dot;

//         USparseMatrix mass_mc = mass_m;
//         USparseMatrix mass_lumped_mc = mass_lumped_m;

//         USparseMatrix mass_fc = mass_f;
//         USparseMatrix mass_lumped_fc = mass_lumped_f;

//         constraint_concentration_mat(rhs_m_c, mass_mc, _constraint_m);
//         constraint_concentration_mat(rhs_m_c, mass_lumped_mc, _constraint_m);

//         constraint_concentration_mat(rhs_f_c, mass_fc, _constraint_f);
//         constraint_concentration_mat(rhs_f_c, mass_lumped_fc, _constraint_f);

//         UVector diag_elem_mc = 1./sum((mass_lumped_mc),1);
//         utopia::USparseMatrix inv_mass_lumped_mc = diag(diag_elem_mc);


//         UVector diag_elem_fc = 1./sum((mass_lumped_fc),1);
//         utopia::USparseMatrix inv_mass_lumped_fc = diag(diag_elem_fc);

//         _S_m_store = const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).getTransferOperator();

//         _S_f_store = const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_2)).getTransferOperator();

//         c_dot_m = utopia::local_zeros(utopia::local_size(c_m));

//         c_dot_f = utopia::local_zeros(utopia::local_size(c_m));

//         USparseMatrix A_stab_mc = *_S_m_store;

//         USparseMatrix A_stab_fc = *_S_f_store;

//         // constraint_concentration_mat(rhs_m_c, A_stab_mc, _constraint_m);
//         // constraint_concentration_mat(rhs_m_c, A_stab_fc, _constraint_m);

//         UVector diag_elem_m = 1./sum((mass_lumped_m),1);
//         utopia::USparseMatrix inv_mass_lumped_m = diag(diag_elem_m);
//         c_dot_m = 1.0 * inv_mass_lumped_m * A_stab_mc * c_m;

//         UVector diag_elem_f = 1./sum((mass_lumped_f),1);
//         utopia::USparseMatrix inv_mass_lumped_f = diag(diag_elem_f);
//         c_dot_f = 1.0 * inv_mass_lumped_f * A_stab_fc * c_f;

//         utopia::UVector _f_m = local_zeros(local_size(c_dot_m));
//         utopia::UVector _f_f = local_zeros(local_size(c_dot_f));




//         //stabilize_coeffiecient(_fe_problem, c_m, c_dot_m, A_stab_mc, mass_m, _f_m);

//         //stabilize_coeffiecient(_fe_problem, c_f, c_dot_f, A_stab_fc, mass_f, _f_f);



//         utopia::UVector r_s_m = rhs_m_t - 1.0 * _f_m;

//         utopia::UVector r_s_f = rhs_f_t - 1.0 * _f_f;


//         // constraint_concentration_vec(rhs_m_c,  r_s_m, _constraint_m);
//         utopia::UVector c_m_tot = c_m + 1.0 * dt * inv_mass_lumped_mc * _f_m;

//         // constraint_concentration_vec(rhs_m_c,  r_s_m, _constraint_m);
//         utopia::disp(_f_m);

//         // constraint_concentration_vec(rhs_m_c,  r_s_m, _constraint_m);
//         utopia::UVector c_f_tot = T * c_m_tot;
 
//         // constraint_concentration_vec(rhs_m_c,  r_s_m, _constraint_m);
    
//         // constraint_concentration_vec(rhs_f_c,  r_s_f, _constraint_f);


       

//         // dlagr_t = utopia::local_zeros(local_size(B).get(0));   

//         // lagr_t = utopia::local_zeros(local_size(dlagr_t));

//         // utopia::UVector rhs_tot = utopia::blocks(r_s_m, r_s_f, dlagr_t);

//         // op->update(_A_t_store);
        
//         // op->apply(rhs_tot, sol_t);

//         // utopia::UVector c_m_tot, c_f_tot;

//         // c_m_tot = local_zeros(local_size(c_dot_m));
//         // c_f_tot = local_zeros(local_size(c_dot_f));

//         // utopia::undo_blocks(sol_t, c_m_tot, c_f_tot, lagr_t);


//         CopyMatrixSolution(c_m);
//         CopyFractureSolution(c_f);

//     }
//     }





// /*    UVector P_minus = sum(transpose(P_m_minus),1);

//     UVector P_max = sum(transpose(P_m_max),1);*/



    
// Real
// FractureAppNConforming::ComputeMaterialProprties(const Elem *elem){
    
//    // _console << "_vector_p.size()"  << _vector_p.size() <<std::endl;

// Real poro=0.0;

//     for(int ll=0; ll<_vector_p.size(); ll++){
//         if (elem->subdomain_id()==_vector_p[ll]) {

//             poro = _vector_value[ll]; 
//             //std::cout<<"poro"<<poro<<std::endl;
//         }
//     }
  
//     return poro;   
// }







  



