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
_stab_matrix(_pp_comm),
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

    }

    ParrotProblem::updateSolution();
    
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
        jacobian.add(1.0,_stab_matrix);
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
        _stab_matrix.vector_mult_add(residual,soln);
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
    
    // initialize stabilization matrix
    _stab_matrix.init(m,n,m_l,n_l,30);
    
    Mat _stab_matrix_petsc=_stab_matrix.mat();
    MatSetOption(_stab_matrix_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    
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

                    _stab_matrix.set(row,col,value);
                    _stab_matrix.add(row,row,maxEntry);
                }
            }
            // else{
            //     _stab_matrix.set(row,col,0.0000);
            // }
        }
        MatRestoreRow(jac_petsc,row,&ncols,&cols,&val);
        MatRestoreRow(jac_tr_petsc,row,&ncols_tr,&cols_tr,&val_tr);
    }
    
    _stab_matrix.close();
        
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
    ParrotProblem::stabilize_coeffiecient();

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
 


    if(this->timeStep()==1){

        int m=dof_map.n_dofs();
        
        int n=dof_map.n_dofs();
        
        int m_l=dof_map.n_local_dofs();
        
        int n_l=dof_map.n_local_dofs();
        
        int nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());

        PetscMatrix<Number> _jac_mat(_pp_comm);

        jmat = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).JacMatrix();

        jmat->init(m,n,m_l,n_l,nnz_x_row);

        jmat->zero();
    
        MatCopy(petsc_mat->mat(), jmat->mat(), DIFFERENT_NONZERO_PATTERN);

        jmat->add(-1.0,_stab_matrix);

         ParrotProblem::determine_dc_bnd_var_id(ParrotProblem::split_string(_dc_var, ' '));

        find_boundary(zero_rows, _dc_boundary_id);
    }

    auto _M =  const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).MassMatrix();
    auto _L =  const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).LumpMassMatrix();
    auto _J = const_cast<StoreOperators&>(this->getUserObjectTempl<StoreOperators>(_userobject_name)).JacMatrix();


    PetscVector<Number> _inv(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _inv.zero();

    PetscVector<Number> _u_dot(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _u_dot.zero();

    PetscVector<Number> _tmp(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _tmp.zero();

    
    NumericVector<Number> &ghosted_solution = *_sys.current_local_solution.get();
    

    //std::vector<numeric_index_type> size_p;
    //size_p.clear();
    // for(int kk=0; kk<dof_map.n_dofs(); kk++){
    //        size_p.push_back(kk);
    // }
  
    std::vector<double>_vec_localize;
    ghosted_solution.localize(_vec_localize);

    _L->get_diagonal(_inv);
    _inv.reciprocal();

    _L->vector_mult(_tmp,ghosted_solution);
    _u_dot.pointwise_mult(_inv,ghosted_solution);

    std::vector<double>_vec_localize_dot;
    // _u_dot.localize(_vec_localize_dot);
    
    int r_start = _M->row_start();
    int r_stop = _M->row_stop();


    PetscMatrix<Number> _J_tr(_pp_comm);   
    _J->get_transpose(_J_tr);

    Mat M_petsc    = _M->mat();
    Mat J_tr_petsc = _J_tr.mat();
    Mat J_petsc    = _J->mat(); 
    Mat D_petsc    = _stab_matrix.mat();



    PetscVector<Number> _m_i(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _m_i.zero();


    PetscVector<Number> ones(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    ones.zero();


    // PetscVector<Number> f(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    // f.zero();

    std::unique_ptr<NumericVector<Number>> f = ghosted_solution.zero_clone();

    ones=1.0;
    
    _M->vector_mult(_m_i,ones);



    // ghosted_solution.print_matlab("ghosted_solution.txt");

    // std::ofstream myfile;

    // std::cout<<"_mesh.processor_id()"<<_mesh.processor_id()<<std::endl;

    // std::string a = "example" + std::to_string(_mesh.processor_id());

    // myfile.open (a);
    // for(int i=0; i<_vec_localize.size(); i++){

    //      //myfile << "Writing this to a file.\n";

    //      myfile << _vec_localize.at(i)<<"\n";

    // }
   
    // myfile.close();

    //exit(1);


    for (int row=r_start; row<r_stop; ++row)
    {
        PetscInt ncols_m, ncols_d, ncols_j_tr, ncols_j; 

        PetscInt const *cols_m; 
        PetscInt const *cols_d; 
        PetscInt const *cols_j_tr; 
        PetscInt const *cols_j; 

        PetscScalar const *val_m; 
        PetscScalar const *val_d; 
        PetscScalar const *val_j_tr;
        PetscScalar const *val_j;
        


        std::vector<Real>f_tmp;
        std::vector<Real>alpha_v;
        


        Real _P_p=0.0; 
        Real _Q_p=0.0; 
        Real _R_p=0.0;
        Real _P_m=0.0; 
        Real _Q_m=0.0; 
        Real _R_m=0.0;



        MatGetRow(D_petsc,row,&ncols_d,&cols_d,&val_d);
        Real max_d = std::abs(val_d[0]);

        for (int p=0; p<ncols_d; p++)
        {
            if(max_d<std::abs(val_d[p]))
                 max_d = std::abs(val_d[p]);
        }

        MatRestoreRow(D_petsc,row,&ncols_d,&cols_d,&val_d);

        MatGetRow(J_tr_petsc,row,&ncols_j_tr,&cols_j_tr,&val_j_tr);
        Real max_j_tr = std::abs(val_j_tr[0]);

        for (int p=0; p<ncols_j_tr; p++)
        {
            if(max_j_tr<std::abs(val_j_tr[p]))
                 max_j_tr = std::abs(val_j_tr[p]);
        }

        MatRestoreRow(J_tr_petsc,row,&ncols_j_tr,&cols_j_tr,&val_j_tr);



        Real f_ij = 0.0;
        Real f_ij_sum = 0.0;
        Real max_p = 0.0;
        Real min_p = 0.0;


        MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m);  
        MatGetRow(J_tr_petsc,row,&ncols_j_tr,&cols_j_tr,&val_j_tr);
        MatGetRow(D_petsc,row,&ncols_d,&cols_d,&val_d);
        MatGetRow(J_petsc,row,&ncols_j,&cols_j,&val_j);

        f_tmp.clear();
        alpha_v.clear();
        //f_tmp.resize(ncols_m);

        
       for (int p=0; p<ncols_m; p++)
        {
            int col=cols_m[p];

            if (row!=col)
            {
                
                Real Dij = val_d[p];

                Real Mij = val_m[p];
                
                Real Kij = -1.0 * val_j[p];
                
                Real Kji = -1.0 * val_j_tr[p];

                if(std::abs(Dij)>max_d) Dij=0.0;

                // if(std::abs(Kji)>max_j_tr) Kji=0.0;
                
                f_ij = Dij * std::max(0.0, Kji) * (ghosted_solution(row) - 1.0 * _vec_localize.at(col));

                //std::cout<<"f_ij_in"<<_vec_localize.at(col)<<std::endl;

                //f_ij = Mij * (_u_dot(row)-_vec_localize_dot.at(col)) - Dij * (ghosted_solution(row) - 1.0 * _vec_localize.at(col));

                f_tmp.push_back(f_ij);
                
                f_ij_sum+=f_ij;

                max_p+=std::max(0.0,f_ij);

                min_p+=std::min(0.0,f_ij);             

            }

            if(p==ncols_m-1) {

                //std::cout<<f_tmp.size()<<ncols_m<<p<<std::endl;

                f_tmp.push_back(f_ij_sum);
            }
        }
  
      
        // for (int ll=0; ll<f_tmp.size(); ll++)
        //     std::cout<<f_tmp.at(ll)<<std::endl;

        

        // _P_p = max_p;
        // _P_m = min_p;
        MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);  
        MatRestoreRow(J_tr_petsc,row,&ncols_j_tr,&cols_j_tr,&val_j_tr);
        MatRestoreRow(D_petsc,row,&ncols_d,&cols_d,&val_d);
        MatRestoreRow(J_petsc,row,&ncols_j,&cols_j,&val_j);

        
        Real u_max = ghosted_solution(row);
        Real u_min = ghosted_solution(row);


        // MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m);
        // for (int p=0; p<ncols_m; ++p){
        //         int col=cols_m[p];
        //         u_max = std::max(u_max,_vec_localize.at(col));
        //     }

        // MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);   


        // MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m);    
        // for (int p=0; p<ncols_m; ++p){
        //         int col=cols_m[p];
        //         u_min = std::min(u_min,_vec_localize.at(col));
        // }

        // MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);

        double v_i = ghosted_solution(row);

        // 
    
        Real dt = static_cast<Transient*>(this->getMooseApp().getExecutioner())->getDT();
        
        //double inv_dt = 1.0/dt;        
        
        //_m_i.scale(inv_dt);

        // Real u_i_max = _m_i(row) * max;
        // Real u_i_min = _m_i(row) * min;

        double max = _m_i(row)/dt * (*std::max_element(std::begin(_vec_localize), end(_vec_localize)) - v_i);
        double min = _m_i(row)/dt * (*std::min_element(std::begin(_vec_localize), end(_vec_localize)) - v_i);   
            
        _Q_p = u_max;
        _Q_m = u_min;

        
        Real p_p_q_p = 0.0;
        if(std::abs(_P_p)>0) p_p_q_p = _Q_p/_P_p;

        Real p_m_q_m = 0.0;
        if(std::abs(_P_m)>0) p_m_q_m = _Q_m/_P_m;

        Real r_i_p = std::max(1.0,p_p_q_p);
        Real r_i_m = std::min(1.0,p_m_q_m);

        auto it = std::find(zero_rows.begin(), zero_rows.end(), row);

        //std::cout<<"zero_rows.size()"<<zero_rows.size()<<std::endl;

        if(it == zero_rows.end()){

            _R_p = std::min(1.0,r_i_p);
            _R_m = std::min(1.0,r_i_m);

        }
        else

        {
            _R_p = 1.0;
            _R_m = 1.0;
            //std::cout<<"ciao"<<std::endl;
        }



        //std::cout<<"max_r"<<max_p<<"min_r"<<min_p<<std::endl;
    
        //Real alpha_ij = 0.0;

        //alpha.clear();
        
       
       MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m);  
        // MatGetRow(J_tr_petsc,row,&ncols_j_tr,&cols_j_tr,&val_j_tr);
        // MatGetRow(J_petsc,row,&ncols_j,&cols_j,&val_j);            
        for (int p=0; p<ncols_m; ++p)
        {
            int col=cols_m[p];            
                
            if (row!=col)
            {
                
                if(f_tmp.at(p)>=0) alpha_v.push_back(_R_p);

                else if (f_tmp.at(p)<0) alpha_v.push_back(_R_m);
                
            }
            if(p==ncols_m-1) 
            {    
                
                if(f_tmp.at(p)>=0) alpha_v.push_back(_R_p);

                else if (f_tmp.at(p)<0) alpha_v.push_back(_R_m);
                
            }           
            
        }
        MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);  
        // MatRestoreRow(J_tr_petsc,row,&ncols_j_tr,&cols_j_tr,&val_j_tr);
        // MatRestoreRow(J_petsc,row,&ncols_j,&cols_j,&val_j);


      
       Real ris = 0;

       for (int ll=0; ll<f_tmp.size(); ll++){
           ris+=f_tmp.at(ll)*alpha_v.at(ll);
       }

       //std::cout<<ris<<std::endl;

       ris *= _inv(row);


       (*f).set(row,ris);


       //ghosted_solution(row)+=ris;

    }

  
    ghosted_solution.add(*f);

   // exit(1);
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

  

