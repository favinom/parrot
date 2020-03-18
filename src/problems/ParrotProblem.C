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
#include "petscmat.h"
#include "chrono"



registerMooseObject("parrotApp", ParrotProblem);

template <>
InputParameters
validParams<ParrotProblem>()
{
    InputParameters params = validParams<FEProblem>();
    params.addRequiredParam<bool>("use_AFC","use_AlgFluxCorr");
    params.addParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    return params;
}


static void getRow(PetscMatrix<Number> & matrix, int const & row, std::vector<Real> & values, std::vector<int> & columns)
{
    Mat const & mat=matrix.mat();
    PetscInt ncol;
    PetscInt const *col;
    PetscScalar const *val;
    MatGetRow(mat,row,&ncol,&col,&val);
    values.resize(ncol);
    columns.resize(ncol);
    for (int i=0; i<ncol; ++i)
    {
        values[i] =val[i];
        columns[i]=col[i];
    }
    MatRestoreRow(mat,row,&ncol,&col,&val);
}


ParrotProblem::ParrotProblem(const InputParameters & parameters) :
FEProblem(parameters),
_pp_comm(_mesh.comm()),
_stab_matrix(_pp_comm),
_res_m(_pp_comm),
_use_afc(getParam<bool>("use_AFC"))
{
 
    if(parameters.isParamValid("operator_userobject"))
    {
        _hasStoreOperatorsUO=true;
        userObjectName=new UserObjectName(getParam<UserObjectName>("operator_userobject"));

  // std::vector<UserObject *> userobjs;
  // theWarehouse().query().condition<AttribSystem>("UserObject").queryInto(userobjs);
  // std::cout<<userobjs.size()<<std::endl;
  // exit(1);
    }
    else
    {
        _hasStoreOperatorsUO=false;
        _storeOperatorsUO=NULL;
    }

    PCCreate(PETSC_COMM_WORLD, &_problem_PC);
    PCSetType(_problem_PC,PCLU);

    _const_jacobian=true;
    _is_stab_matrix_assembled=false;
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

    if (_hasStoreOperatorsUO)
    {
        _storeOperatorsUO=&getUserObject<StoreOperators>(userObjectName[0]);
        PetscMatrix<Number> * & localmatrix=_storeOperatorsUO[0].StabMatrix();
        localmatrix=&_stab_matrix;
    }


    NonlinearSystemBase & _nl = this->getNonlinearSystemBase();
    const DofMap & dof_map = _nl.dofMap(); 
    _sol_vec                = _storeOperatorsUO->SolVec();
    _sol_vec->init(dof_map.n_dofs(), dof_map.n_local_dofs());

    std::cout<<"END ParrotProblem::initialSetup"<<std::endl;
    
};


void
ParrotProblem::solve()
{


    std::cout<<"I am in FEProblemBase::solve()"<<std::endl;

    auto t_start = std::chrono::high_resolution_clock::now();


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
    _fail_next_linear_convergence_check = false;

    if (_solve) _nl->solve();
    

    if (_solve)
    {

      if(this->timeStep()==1)_nl->update();
      
      else update_sol();
       
    }

}

void ParrotProblem::timestepSetup()
{
    std::cout<<"BEGIN ParrotProblem::timestepSetup"<<std::endl;
    
    FEProblem::timestepSetup();
    
    PetscErrorCode ierr;
    Moose::PetscSupport::petscSetOptions(*this);

    _sol_vec->zero();
    
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
    _ksp_ptr->userObjectName=userObjectName;
    _ksp_ptr->fe_problem = this;
    
    std::cout<<"END ParrotProblem::timestepSetup"<<std::endl;
};

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
        

    if(this->timeStep()==1){

        NonlinearImplicitSystem * a;

        NonlinearSystemBase & nl_sys = this->getNonlinearSystemBase();

        DofMap const & dof_map = nl_sys.dofMap();

        FEProblemBase::computeResidualSys(a[0],soln,residual);

        _res_m.init(dof_map.n_dofs(), dof_map.n_local_dofs());

        if (_use_afc && _is_stab_matrix_assembled)
        {
        _stab_matrix.vector_mult_add(residual,soln);

        //residual.print_matlab("residual_moose.m");
        }

        _poro_lumped = _storeOperatorsUO->PoroLumpMassMatrix();

        _dirichlet_bc = _storeOperatorsUO->BcVec();

        _value_dirichlet_bc = _storeOperatorsUO->ValueBcVec();

    }
    else
    {
     
        auto t_start = std::chrono::high_resolution_clock::now();

        std::cout<<"start computing?\n";

        _res_m.zero();
    
        _poro_lumped->vector_mult(_res_m, soln);

        residual.pointwise_mult(_res_m,*_dirichlet_bc);

        Real inv_dt = 1.0/this->dt();
        
        residual.scale(inv_dt);

        residual.add(*_value_dirichlet_bc);

        auto t_end = std::chrono::high_resolution_clock::now();

        std::cout<<"done computing?\n";
        std::cout<<"computing time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";

    }

 
    
}


void
ParrotProblem::computeStabilizationMatrix(SparseMatrix<Number> & jacobian)
{
    std::cout<<"START ParrotProblem::computeStabilizationMatrix\n";
    PetscMatrix<Number> & jac_PM=dynamic_cast<PetscMatrix<Number> &> (jacobian);
   // Declare a PetscMatrix that will contain the transpose
    PetscMatrix<Number> jac_tr_PM(_pp_comm);
    // Transpose of the Jacobian
    jacobian.get_transpose (jac_tr_PM);
        
    int rb=jac_PM.row_start();
    int re=jac_PM.row_stop();

    std::cout<<"get dof map\n";
    DofMap const & dof_map = this->getNonlinearSystemBase().dofMap();
    std::cout<<"done\n";

    std::cout<<"attach dofmap\n";
    _stab_matrix.attach_dof_map(dof_map);
    std::cout<<"done\n";

    std::cout<<"init matrix\n";
    _stab_matrix.init();//m,n,m_l,n_l,30);
    std::cout<<"done\n";

    std::cout<<"zeroing\n";
    _stab_matrix.zero();
    std::cout<<"done\n";
    //MatSetOption(_stab_matrix_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    
    std::vector<Real> values;
    std::vector<Real> values_tr;
    std::vector<int> columns;
    std::vector<int> columns_tr;

    std::cout<<"looping\n";
    for (int row=rb; row<re; ++row)
    {
        getRow(jac_PM,row,values,columns);
        getRow(jac_tr_PM,row,values_tr,columns_tr);
        
        if (columns.size()!=columns_tr.size())
        {
            std::cout<<"ncols!=ncols_tr\n";
            exit(1);
        }
        
        for (int p=0; p<columns.size(); ++p)
        {
            if (columns[p]!=columns_tr[p])
            {
                std::cout<<"cols[p]!=cols_tr[p]"<<std::endl;
                exit(1);
            }
            
            int col=columns[p];
            if (row!=col)
            {
                // we are in a extradiagonal
                Real Aij=values[p];
                Real Aji=values_tr[p];
                
                Real maxEntry=std::max(Aij,Aji);
                
                if (maxEntry>-1e-15)
                {
                    Real value=-1.0*maxEntry;

                    _stab_matrix.set(row,col,value);
                    _stab_matrix.add(row,row,maxEntry);
                }
            }
        }
    }
    std::cout<<"done\n";

    std::cout<<"closing\n";
     _stab_matrix.close();
    std::cout<<"done\n";

    _is_stab_matrix_assembled=true;

    std::cout<<"STOP ParrotProblem::computeStabilizationMatrix\n";
}


void
ParrotProblem::update_sol()
{
    std::cout<<"ParrotProblem::update_sol"<<std::endl;

    MooseVariableFEBase  & main_var = this->getVariable(0, "CM", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

    System & main_sys = main_var.sys().system();

    NumericVector<Number> * main_solution = main_sys.solution.get();


    _sol_vec                = _storeOperatorsUO->SolVec();

    //_sol_vec->print_matlab("sol_moose.m");

    { 
        
        for (const auto & node : this->es().get_mesh().local_node_ptr_range())

        {
            for (unsigned int comp = 0; comp < node->n_comp(main_sys.number(), main_var.number()); comp++)

            {

                //const dof_id_type proj_index = node->dof_number(_nl.number(), sol_var.number(), comp);

                const dof_id_type to_index = node->dof_number(main_sys.number(), main_var.number(), comp);

                //const dof_id_type to_index_c = node->dof_number(aux_sys.number(), aux_var_c.number(), comp);


                main_solution->set(to_index, (*_sol_vec)(to_index));

                //aux_solution->set(to_index_c, correction(proj_index));
            }

        }
    }



    main_solution->close();
    //aux_solution->close();
    main_sys.update();



 
    
  

}


