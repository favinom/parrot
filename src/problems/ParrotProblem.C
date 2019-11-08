//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ParrotProblem.h"

#include "libmesh/petsc_nonlinear_solver.h"

#include "petscmat.h"

#include "Transient.h"

#include <iostream>

#include "TransientMultiApp.h"

registerMooseObject("parrotApp", ParrotProblem);

template <>
InputParameters
validParams<ParrotProblem>()
{
    InputParameters params = validParams<FEProblem>();
    params.addRequiredParam<bool>("use_AFC","use_AlgFluxCorr");
    params.addRequiredParam<bool>("change_sol","change_sol");
    params.addParam<std::string>("dc_boundaries", "-1", "Dirichlet Boundary ID");
    params.addParam<std::string>("dc_variables" , "-1", "Variable to which given BC_id applies");
    //params.addRequiredParam<UserObjectName>("operator_userobject_problem","The userobject that stores our operators");
    return params;
}

ParrotProblem::ParrotProblem(const InputParameters & parameters) :
FEProblem(parameters),
_pp_comm(_mesh.comm()),
// _stab_matrix(_pp_comm),
_use_afc(getParam<bool>("use_AFC")),
_change_sol(getParam<bool>("change_sol")),
_dc_var(getParam<std::string>("dc_variables"))
//_operator_storage(getUserObject<StoreOperators>("operator_userobject_problem"))
{
    _const_jacobian=true;
    PCCreate(PETSC_COMM_WORLD, &_problem_PC);
    PCSetType(_problem_PC,PCLU);
    
    _is_stab_matrix_assembled=false;

    std::vector<std::string> tmp = split_string(parameters.get<std::string>("dc_boundaries"), ' ');
    for(auto str_tmp=tmp.begin(); str_tmp != tmp.end(); str_tmp++)
    {
        _dc_boundary_id.push_back(atoi(str_tmp->c_str()));
    }

   

}

void ParrotProblem::initialSetup()
{
    std::cout<<"BEGIN ParrotProblem::initialSetup"<<std::endl;
    
    FEProblem::initialSetup();
    
    Moose::PetscSupport::petscSetOptions(*this);
    
    PetscErrorCode ierr;
    
    PetscNonlinearSolver<Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<Number> *>((*_nl_sys).nonlinearSolver());
    
    SNES snes = petsc_solver->snes();
    KSP ksp;
    SNESType ttype;
    KSPType type;
    ierr = SNESGetKSP(snes,&ksp);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);
    std::cout<< "SNES type from PARROTPROBLEM exec:"<<ttype<<std::endl;
    ierr = KSPGetType(ksp, &type);
    std::cout<< "KSP  type from PARROTPROBLEM exec:"<<type<<std::endl;
    _factorized = 0 ;

    if(!this->hasUserObject(_userobject_name)){

    _console << "Add UserObject" << std::endl;

    std::string class_name = "StoreOperators";
    auto params = this->getMooseApp().getFactory().getValidParams(class_name);
    params.set<bool>("use_displaced_mesh") = false;
    params.set<ExecFlagEnum>("execute_on") = "initial";

    this->addUserObject("StoreOperators", _userobject_name, params);
    }
    
    std::cout<<"END ParrotProblem::initialSetup"<<std::endl;
    
};

void ParrotProblem::timestepSetup()
{
    std::cout<<"BEGIN ParrotProblem::timestepSetup"<<std::endl;
    
    FEProblem::timestepSetup();
    
    PetscErrorCode ierr;
    Moose::PetscSupport::petscSetOptions(*this);
    
    PetscNonlinearSolver<Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<Number> *>((*_nl_sys).nonlinearSolver());
    SNES snes = petsc_solver->snes();
    KSP ksp;
    PC pc;
    SNESType ttype;
    KSPType type;
    PCType ptype;
    
    ierr = SNESGetKSP(snes,&ksp);
    ierr = KSPGetPC(ksp,&pc);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);
    ierr = KSPGetType(ksp, &type);
    ierr = PCGetType(pc, &ptype);
    
    std::cout<< "SNES type da passo steady exec:"<<ttype<<std::endl;
    std::cout<< "KSP type da passo steady exec:"<<type<<std::endl;
    std::cout<< "PC type da passo steady exec:"<<ptype<<std::endl;
    
    _ksp_ptr = (KSP_PARROT *)ksp->data;
    (_ksp_ptr[0].local_pc)=&_problem_PC;
    (*_ksp_ptr).factorized=&_factorized;
    
    std::cout<<"END ParrotProblem::timestepSetup"<<std::endl;
};

void ParrotProblem::solve()
{
    std::cout<<"BEGIN ParrotProblem::solve"<<std::endl;
    
    Moose::perf_log.push("solve()", "Execution");
    
#ifdef LIBMESH_HAVE_PETSC
    Moose::PetscSupport::petscSetOptions(*this); // Make sure the PETSc options are setup for this app
#endif
    
    Moose::setSolverDefaults(*this);
    
    // Setup the output system for printing linear/nonlinear iteration information
    initPetscOutput();
    
    possiblyRebuildGeomSearchPatches();
    
    // reset flag so that linear solver does not use
    // the old converged reason "DIVERGED_NANORINF", when
    // we throw  an exception and stop solve
    //_fail_next_linear_convergence_check = false;
    
    if (_solve)
    {

        _nl->solve();
  
    }
    
    if (_solve)
    {
        _nl->update();

        ParrotProblem::updateSolution();

    }


    
    // sync solutions in displaced problem
    //if (_displaced_problem)
    //    _displaced_problem->syncSolutions();
    
    Moose::perf_log.pop("solve()", "Execution");
    
    std::cout<<"END ParrotProblem::solve"<<std::endl;
}

void
ParrotProblem::computeJacobianSys(NonlinearImplicitSystem & /*sys*/,
                                  const NumericVector<Number> & soln,
                                  SparseMatrix<Number> & jacobian)
{
    FEProblemBase::computeJacobian(soln, jacobian);
    
    if (_use_afc==true && _is_stab_matrix_assembled==false)
    {
        computeStabilizationMatrix(jacobian);
        jacobian.add(1.0,*_stab_matrix);
    }
}

void
ParrotProblem::computeResidualSys(NonlinearImplicitSystem & /*sys*/,
                                  const NumericVector<Number> & soln,
                                  NumericVector<Number> & residual)
{
    
    //std::cout<<"CIAO sono in Residuo"<<std::endl;
    
    NonlinearImplicitSystem * a;
    FEProblemBase::computeResidualSys(a[0],soln,residual);
    
    PetscVector<Number> ciao(_pp_comm);
    
    if (_use_afc && _is_stab_matrix_assembled)
    {
        //_stab_matrix.vector_mult(ciao,soln);
        _stab_matrix->vector_mult_add(residual,soln);
    }
    
    //std::cout<<"ho finito Residuo"<<std::endl;
    
}


void
ParrotProblem::computeStabilizationMatrix(SparseMatrix<Number> & jacobian)
{
    // Cast pf the jacobian to a Petsc matrix
    PetscMatrix<Number> & jac_PM=dynamic_cast<PetscMatrix<Number> &> (jacobian);


    
    // Declare a PetscMatrix that will contain the transpose
    PetscMatrix<Number> jac_tr_PM(_pp_comm);
    // Transpose of the Jacobian
    jacobian.get_transpose (jac_tr_PM);
    
    // Get the petsc matrix (Mat) to call MatGetRow,
    //the only function missing in the wrapper
    Mat jac_petsc=jac_PM.mat();
    Mat jac_tr_petsc=jac_tr_PM.mat();
    
    int m=jac_PM.m();
    int n=jac_PM.n();
    
    int m_l=jac_PM.local_m();
    int n_l=jac_PM.local_n();
    
    int rb=jac_PM.row_start();
    int re=jac_PM.row_stop();

    _stab_matrix = const_cast<StoreOperators&>(getUserObjectTempl<StoreOperators>(_userobject_name)).StabMatrix();
    
    // initialize stabilization matrix
    (*_stab_matrix).init(m,n,m_l,n_l,30);
    
    //Mat _stab_matrix_petsc=_stab_matrix.mat();
    //MatSetOption(_stab_matrix_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    
    MatCopy(jac_PM.mat(),  (*_stab_matrix).mat(), DIFFERENT_NONZERO_PATTERN);
    (*_stab_matrix).zero();

    for (int row=rb; row<re; ++row)
    {
        PetscInt ncols;
        PetscInt const *cols;
        PetscScalar const *val;
        MatGetRow(jac_petsc,row,&ncols,&cols,&val);
        PetscInt ncols_tr;
        PetscInt const *cols_tr;
        PetscScalar const *val_tr;
        MatGetRow(jac_tr_petsc,row,&ncols_tr,&cols_tr,&val_tr);
        
        if (ncols!=ncols_tr)
        {
            std::cout<<"ncols!=ncols_tr\n";
            exit(1);
        }
        
        for (int p=0; p<ncols; ++p)
        {
            if (cols[p]!=cols_tr[p])
            {
                std::cout<<"cols[p]!=cols_tr[p]"<<std::endl;
                exit(1);
            }
            
            int col=cols[p];
            

            
            if (row!=col)
            {
                // we are in a extradiagonal
                Real Aij=val[p];
                Real Aji=val_tr[p];
                
                Real maxEntry=std::max(Aij,Aji);
                
                if (maxEntry>-1e-15)
                {
                    Real value=-1.0*maxEntry;

                    (*_stab_matrix).set(row,col,value);
                    (*_stab_matrix).add(row,row,maxEntry);
                }
            }
            // else{
            //     _stab_matrix.set(row,col,0.0000);
            // }
        }
        MatRestoreRow(jac_petsc,row,&ncols,&cols,&val);
        MatRestoreRow(jac_tr_petsc,row,&ncols_tr,&cols_tr,&val_tr);
    }
    
     (*_stab_matrix).close();
        
    _is_stab_matrix_assembled=true;

    // ParrotProblem::stabilize_coeffiecient();
    
}

bool
 ParrotProblem::shouldUpdateSolution()
 {
    if (_change_sol) {
         //ParrotProblem::stabilize_coeffiecient();
         return true;
    }
    else 
        return false;
 }


bool
 ParrotProblem::updateSolution()
 {
    //ParrotProblem::stabilize_coeffiecient();

    return true;
 }


void
ParrotProblem::stabilize_coeffiecient(/*NumericVector<Number> & vec_solution,
                                      NumericVector<Number> & ghosted_solution*/){



    auto &_sys = this->es().get_system<TransientNonlinearImplicitSystem>("nl0");

    //_stab_matrix.print_matlab("original.txt");

    NonlinearSystemBase & _nl = this->getNonlinearSystemBase();
    
    PetscMatrix<Number> *petsc_mat = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_sys.matrix);

    //petsc_mat->print_matlab("petsc_mat.txt");
    

    const DofMap & dof_map = this->getNonlinearSystemBase().dofMap();


    int m=dof_map.n_dofs();
    
    int n=dof_map.n_dofs();
    
    int m_l=dof_map.n_local_dofs();
    
    int n_l=dof_map.n_local_dofs();
        

    int nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());

    if(this->timeStep()==1){


        PetscMatrix<Number> _jac_mat(_pp_comm);

        jmat = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).JacMatrix();

        jmat->init(m,n,m_l,n_l,nnz_x_row);

        
    
        MatCopy(petsc_mat->mat(), jmat->mat(), DIFFERENT_NONZERO_PATTERN);

        jmat->zero();

        jmat->add(-1.0,*_stab_matrix);

        ParrotProblem::determine_dc_bnd_var_id(ParrotProblem::split_string(_dc_var, ' '));

        find_boundary(zero_rows, _dc_boundary_id);
    }

    auto _M  = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).MassMatrix();
    auto _L  = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).LumpMassMatrix();
    auto _J  = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).JacMatrix();
    auto _PL = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).PoroLumpMassMatrix();
    auto _PM = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).PoroMassMatrix();


    PetscVector<Number> _inv(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _inv.zero();

    PetscVector<Number> _inv_p(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _inv_p.zero();

    // PetscVector<Number> _u_dot(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    // _u_dot.zero();

    PetscVector<Number> _tmp(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _tmp.zero();

    
    NumericVector<Number> &ghosted_solution = *_sys.current_local_solution.get();
    
  
    std::vector<double>_vec_localize;
    ghosted_solution.localize(_vec_localize);

    _L->get_diagonal(_inv);
    _inv.reciprocal();

    _PL->get_diagonal(_inv_p);
    _inv_p.reciprocal();

    NumericVector<Number> * _u_dot_moose = _nl.solutionUDot();

    PetscVector<Number> &_u_dot = dynamic_cast<libMesh::PetscVector<libMesh::Number>& >(*_u_dot_moose);

    std::vector<double>_vec_localize_dot;
    _u_dot.localize(_vec_localize_dot);


    NumericVector<Number> * local_vector;
    std::unique_ptr<NumericVector<Number>> local_vector_built;
    local_vector_built = NumericVector<Number>::build(dof_map.comm());
    local_vector = local_vector_built.get();
    local_vector->init(dof_map.n_dofs(), false, SERIAL);
    ghosted_solution.localize(*local_vector,dof_map.get_send_list());
    local_vector->close();

    local_vector->print_matlab("local_vector.m");
    
    int r_start = _PM->row_start();
    int r_stop  = _PM->row_stop();


    PetscMatrix<Number> _J_tr(_pp_comm);   
    _J->get_transpose(_J_tr);

    Mat PM_petsc    = _PM->mat();
    Mat L_petsc     = _L->mat();
    Mat D_petsc     = (*_stab_matrix).mat();


    PetscVector<Number> _R_p(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _R_p.zero();

    PetscVector<Number> _R_m(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _R_m.zero();

    PetscVector<Number> _a_bar(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _a_bar.zero();

    PetscVector<Number> _m_i(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _m_i.zero();

    PetscVector<Number> _ones(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _ones.zero();


    _ones.add(1.0);


    PetscMatrix<Number> _f_mat(_pp_comm);
    _f_mat.init(m,n,m_l,n_l,nnz_x_row);     
    Mat _f_mat_petsc =_f_mat.mat();
    MatCopy(petsc_mat->mat(), _f_mat_petsc, DIFFERENT_NONZERO_PATTERN);
    _f_mat.zero();
    //MatSetOption(_f_mat_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   

    PetscMatrix<Number> _alpha_mat(_pp_comm);
    _alpha_mat.init(m,n,m_l,n_l,nnz_x_row);
    Mat _alpha_mat_petsc = _alpha_mat.mat();
    _alpha_mat.zero();
    MatSetOption(_alpha_mat_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    

    std::unique_ptr<NumericVector<Number>> f = ghosted_solution.zero_clone();

    
    _L->vector_mult(_m_i,_ones);


    Real dt = static_cast<Transient*>(this->getMooseApp().getExecutioner())->getDT();



    for (int row=r_start; row<r_stop; ++row)
    {
        PetscInt ncols_m, ncols_d, ncols_l; 

        PetscInt const *cols_m; 
        PetscInt const *cols_d; 
        PetscInt const *cols_l; 


        PetscScalar const *val_m; 
        PetscScalar const *val_d; 
        PetscScalar const *val_l;
        



        Real _P_p=0.0; 
        Real _Q_p=0.0; 

        Real _P_m=0.0; 
        Real _Q_m=0.0; 

        

        Real f_ij = 0.0;
        Real f_ij_sum = 0.0;



        MatGetRow(PM_petsc,row,&ncols_m,&cols_m,&val_m);  
        MatGetRow(D_petsc,row,&ncols_d,&cols_d,&val_d);


        
        for (int p=0; p<ncols_m; p++)
        {
            int col=cols_m[p];

            if (row!=col)
            {
                
                Real Dij  = val_d[p];

                Real PMij = val_m[p];
                
                auto v_i = ghosted_solution(row);

                auto v_d_i = _u_dot(row);

                if(_vec_localize.at(col)>0.010000001){
                    std::cout<<"col_p"<<col<<" and value_p "<<_vec_localize.at(col)<<std::endl;
                }

                if(_vec_localize.at(col)<-1.0e-5){
                    std::cout<<"col_m"<<col<<" and value_m "<<_vec_localize.at(col)<<std::endl;
                }

                if(_vec_localize_dot.at(col)>0.010000001){
                    std::cout<<"col_p"<<col<<" and value_p "<<_vec_localize.at(col)<<std::endl;
                }

                if(_vec_localize_dot.at(col)<-1.0e-5){
                    std::cout<<"col_m"<<col<<" and value_m "<<_vec_localize.at(col)<<std::endl;
                }

                f_ij = PMij * (v_d_i  -_vec_localize_dot.at(col)) - 1.0 * Dij * (v_i - _vec_localize.at(col));


                double check = (v_i - _vec_localize.at(col)) * f_ij;

                if(check==0) _f_mat.set(row, col,0.0);
                else _f_mat.set(row, col, f_ij);

                _P_p+=std::max(0.0,f_ij);

                _P_m+=std::min(0.0,f_ij);             

            }

        }



        MatRestoreRow(PM_petsc,row,&ncols_m,&cols_m,&val_m);  
        MatRestoreRow(D_petsc,row,&ncols_d,&cols_d,&val_d);


        
        Real u_max = ghosted_solution(row);

        Real u_min = ghosted_solution(row);

        std::vector<Real> values;
        values.clear();


        MatGetRow(PM_petsc,row,&ncols_m,&cols_m,&val_m);

        for (int p=0; p<ncols_m; ++p){
                int col=cols_m[p];
                values.push_back(_vec_localize.at(col));
            }

        MatRestoreRow(PM_petsc,row,&ncols_m,&cols_m,&val_m);   

        double v_i = ghosted_solution(row);

        double m_i = _m_i(row);

        _Q_p = m_i/dt * (*std::max_element(std::begin(values), end(values)) - v_i);

        _Q_m = m_i/dt * (*std::min_element(std::begin(values), end(values)) - v_i);   
            
        
        Real value_p = 0.0;

        if(std::abs(_P_p)>0.0) {

            value_p = _Q_p/_P_p;
        }

        Real value_m = 0.0;

        if(std::abs(_P_m)>0.0) {

            value_m = _Q_m/_P_m;
        }

        Real entry_p = std::min(1.0,value_p);

        Real entry_m = std::min(1.0,value_m);

        auto it = std::find(zero_rows.begin(), zero_rows.end(), row);

        //std::cout<<"zero_rows.size()"<<zero_rows.size()<<std::endl;

        if(it == zero_rows.end()){

            _R_p.set(row,entry_p);

            _R_m.set(row,entry_m);
         

        }
        else

        {
            _R_p.set(row,1.0);

            _R_m.set(row,1.0);
            //std::cout<<"ciao"<<std::endl;
        }
    }

     
    _f_mat.close();
   
    _R_m.close();

    _R_p.close();

    _f_mat.print_matlab("f_d.m");

    std::vector<double>_vec_localize_r_p;
    _R_p.localize(_vec_localize_r_p);

    if (*std::max_element(std::begin(_vec_localize_r_p), end(_vec_localize_r_p))>1.0) {
        std::cout<<"max_p="<<*std::max_element(std::begin(_vec_localize_r_p), end(_vec_localize_r_p))<<std::endl;
        exit(1);
    }

        if (*std::min_element(std::begin(_vec_localize_r_p), end(_vec_localize_r_p))<0.0) {
        std::cout<<"minx_p="<<*std::min_element(std::begin(_vec_localize_r_p), end(_vec_localize_r_p))<<std::endl;
        exit(1);
    }


    std::vector<double>_vec_localize_r_m;
    _R_m.localize(_vec_localize_r_m);


    if (*std::max_element(std::begin(_vec_localize_r_m), end(_vec_localize_r_m))>1.0) {
        std::cout<<"max_m="<<*std::max_element(std::begin(_vec_localize_r_m), end(_vec_localize_r_m))<<std::endl;
        exit(1);
    }


    if (*std::min_element(std::begin(_vec_localize_r_m), end(_vec_localize_r_m))<0.0) {
        std::cout<<"minx_m="<<*std::min_element(std::begin(_vec_localize_r_m), end(_vec_localize_r_m))<<std::endl;
        exit(1);
    }




    _R_m.print_matlab("R_m.m");

    _R_p.print_matlab("R_p.m");



    for (int row=r_start; row<r_stop; ++row)
    {
        PetscInt ncols_m, ncols_f, ncols_l; 

        PetscInt const *cols_m, *cols_f, *cols_l; 
     
        PetscScalar const *val_m, *val_f, *val_l; 

        std::vector<double> values;
        values.clear();
        

        MatGetRow(PM_petsc,row,&ncols_m,&cols_m,&val_m); 
        for (int p=0; p<ncols_m; ++p){
                int col=cols_m[p];
                values.push_back(_vec_localize.at(col));
            }
        MatRestoreRow(PM_petsc,row,&ncols_m,&cols_m,&val_m);   

        
        MatGetRow(_f_mat.mat(),row,&ncols_f,&cols_f,&val_f);  
        MatGetRow(PM_petsc,row,&ncols_m,&cols_m,&val_m); 

        double v_i = ghosted_solution(row);    

        for (int p=0; p<ncols_m; ++p)
        {
            int col=cols_m[p];

            if(row!=col) {

                Real ris = 0.0;

                if(val_f[p]>0.0) {
                
                    double entry_p = std::min(_R_p(row),_vec_localize_r_m.at(col));

                    ris = val_f[p] * entry_p;


                    if(std::abs(ris)>1.0){
                        std::cout<<"c_p"<<col<<" and value_p "<<ris<<"  and "<<_R_p(row)<<" and "<<val_f[p]<<std::endl;
                    }


                    // if(std::abs(*std::max_element(std::begin(values), end(values)) - v_i)<1.0e-7){
                    //     _alpha_mat.set(row, col, 0.0);
                    // }
                    // else{
                        _alpha_mat.set(row, col, ris);
                    //}
                 
                    
                }

                else if(val_f[p]<0.0) 
                {

                    double entry_m = std::min(_R_m(row),_vec_localize_r_p.at(col));

                    ris = val_f[p] * entry_m;

                    // if(std::abs(*std::min_element(std::begin(values), end(values)) - v_i)<1.0e-7){

                    //     _alpha_mat.set(row, col, 0.0);

                    // }
                    // else{


                    if(std::abs(ris)>1.0){
                        std::cout<<"c_m"<<col<<" and value_m "<<ris<<" and "<<_R_m(row)<<" and "<<val_f[p]<<std::endl;
                    }


                        _alpha_mat.set(row, col, ris);
                    
                    //}
                    
                }
                
                else 
                {

                    _alpha_mat.set(row, col, 0.0);
                  
                }
            }
                
        }

        MatRestoreRow(PM_petsc,row,&ncols_m,&cols_m,&val_m);  
 
        MatRestoreRow(_f_mat.mat(),row,&ncols_f,&cols_f,&val_f); 


    }


    _alpha_mat.close();


    _alpha_mat.print_matlab("alpha.m");

    _alpha_mat.vector_mult(_a_bar,_ones);


    for (int row=r_start; row<r_stop; ++row)
    {
          
        auto it = std::find(zero_rows.begin(), zero_rows.end(), row);

        if(it == zero_rows.end()){

            auto f_bar = _a_bar(row) * _inv(row) * dt;

            (*f).set(row,f_bar);
        }
        else{

           (*f).set(row,0.0);
        }


    }


             
        
    (*f).close();

    // (*f).print_matlab("f.m");

    ghosted_solution.add(*f);

    // ghosted_solution.print_matlab("sol.m");

    PetscVector<Number> &f_p = dynamic_cast<libMesh::PetscVector<libMesh::Number>& >(ghosted_solution);
    set_solution(f_p);

}

void
ParrotProblem::find_boundary(std::vector<int> &zero_rows, std::vector<int> &_dc_boundary_id){

  ConstBndNodeRange & bnd_nodes = *_mesh.getBoundaryNodeRange();
    unsigned int i = 0;
    // std::cout<<"_dc_variables_id"<< _dc_variables_id[0].size()<<std::endl;
    // std::cout<<"_dc_boundary_id"<< _dc_boundary_id.size()<<std::endl;
    for(auto boundary = _dc_boundary_id.begin(); boundary != _dc_boundary_id.end(); ++boundary, i++)
      {
        // iterate just over boundary nodes
            for (const auto & bnode : bnd_nodes)
            {
                  libMesh::Node * current_node = bnode->_node;

                  // check if node is in active boundary list
                  if (_mesh.isBoundaryNode(current_node->id(), *boundary))
                  {
                    // loop over all variables at this node

                    for (auto v = 0; v < this->getNonlinearSystemBase().nVariables(); v++)
                    {
                      const Variable & var = this->getNonlinearSystem().sys().variable(v);
                      unsigned int var_num = var.number();
                        //std::cout<<"nnnnnnn"<< var_num <<std::endl;

                      // see if this variable has any dofs at this node
                      if (current_node->n_dofs(this->getNonlinearSystem().number(), var_num) > 0)
                      {
                        // check if given variable has BC on node

                        if(std::find(_dc_variables_id[i].begin(), _dc_variables_id[i].end(), var_num) != _dc_variables_id[i].end())
                        {

                          // different components are not supported by moose at the moment...
                          //std::cout<<"kkkkkkkk"<< std::endl;
                          zero_rows.push_back(
                              current_node->dof_number(this->getNonlinearSystem().number(), var_num, 0));
                        }
                    }
                }
            } 
        }
    }

    //std::cout<<"zero_rows"<< zero_rows.size()<<std::endl;
}


void 
ParrotProblem::determine_dc_bnd_var_id(const std::vector<std::string> & BC_var){
    // automatic fill-in
     NonlinearSystemBase & _nl = this->getNonlinearSystemBase();

    std::vector<int> vec(_nl.nVariables());

    std::iota(vec.begin(), vec.end(), 0);

    unsigned int i;

    auto str_tmp = BC_var.begin();

    PetscFunctionBegin;
    // going over all BC_ids
    for(i = 0; str_tmp != BC_var.end(); i++, str_tmp++)
    {
        std::vector<std::string> tmp = ParrotProblem::split_string(*str_tmp, '-');

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
    std::cout<<" ------ BC CONDITIONS  ------ \n";
    unsigned int t = 0;
    //std::cout<<"_dc_variables_id.begin()"<<_dc_variables_id.size()<<std::endl;
    for(auto i = _dc_variables_id.begin(); i != _dc_variables_id.end();  t++, i++)
    {
        std::cout<<"\n BC_id:  "<< _dc_boundary_id[t] << "   var_ids:  ";
        std::for_each(i->begin(), i->end(), [](int i){ std::cout << i << "  " ; });
    }

}


    

     
std::vector<std::string>
ParrotProblem::split_string(const std::string & s, char delim)
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

void ParrotProblem::set_solution(PetscVector<Number> &correction)
{
    // copy projected solution into target es

   


    
    MooseVariableFEBase  & aux_var = getVariable(0, "correction", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

    MooseVariableFEBase  & sol_var = getVariable(0, "CM", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
    
    // solution of the original system
    System & aux_sys = aux_var.sys().system();



    NumericVector<Number> * aux_solution = aux_sys.solution.get();


    NonlinearSystemBase & _nl = this->getNonlinearSystemBase();

    //NumericVector<Number> & from_solution = *ls.solution;
    
  
    //LinearImplicitSystem & ls = _es.get_system<LinearImplicitSystem>("Diffusion");

  
    { // loop through our local elements and set the solution from the projection
        
        for (const auto & node : es().get_mesh().local_node_ptr_range())

        {
            for (unsigned int comp = 0; comp < node->n_comp(aux_sys.number(), aux_var.number()); comp++)

            {

            //std::cout<<"uno"<<std::endl;

                const dof_id_type proj_index = node->dof_number(_nl.number(), sol_var.number(), comp);
            
            //std::cout<<"due"<<std::endl;

                const dof_id_type to_index = node->dof_number(aux_sys.number(), aux_var.number(), comp);

            //std::cout<<"tre"<<std::endl;

             //

                aux_solution->set(to_index, correction(proj_index));
            }

        }
    }

  //auto &aux_sys = _fe_problem.getAuxiliarySystem();




  aux_solution->close();
  aux_sys.update();

  //ExodusII_IO (_fe_problem.es().get_mesh()).write_equation_systems("matrix_c.e", _fe_problem.es());

    
}
  


    // for (int row=r_start; row<r_stop; ++row)
    // {
    //     PetscInt ncols_a, ncols_f; 

    //     PetscInt const *cols_a, *cols_f; 
     
    //     PetscScalar const *val_a, *val_f; 

          
    //     MatGetRow(_alpha_mat.mat(),row,&ncols_a,&cols_a,&val_a); 

    //     MatGetRow(_f_mat.mat(),row,&ncols_f,&cols_f,&val_f);  

    //     double value=0.0;

    //     for (int p=0; p<ncols_a; ++p)
    //     {
    //         int col=cols_a[p];

    //         if(row!=col) {

    //             value+=val_a[p]*val_f[p];
    //         }
    //     }

    //     auto it = std::find(zero_rows.begin(), zero_rows.end(), row);

    //     if(it == zero_rows.end()){

    //         auto f_bar = value * _inv_p(row);

    //         (*f).set(row,f_bar);
    //     }
    //     else{
    //        (*f).set(row,0.0);
    //     }

    //     MatRestoreRow(_alpha_mat.mat(),row,&ncols_a,&cols_a,&val_a);  
    //     MatRestoreRow(_f_mat.mat(),row,&ncols_f,&cols_f,&val_f); 
    // }
