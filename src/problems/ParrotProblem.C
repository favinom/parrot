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

registerMooseObject("parrotApp", ParrotProblem);

template <>
InputParameters
validParams<ParrotProblem>()
{
    InputParameters params = validParams<FEProblem>();
    params.addRequiredParam<bool>("use_AFC","use_AlgFluxCorr");
    return params;
}

ParrotProblem::ParrotProblem(const InputParameters & parameters) :
FEProblem(parameters),
_pp_comm(_mesh.comm()),
_stab_matrix(_pp_comm),
_use_afc(getParam<bool>("use_AFC"))
{
    _const_jacobian=true;
    PCCreate(PETSC_COMM_WORLD, &_problem_PC);
    PCSetType(_problem_PC,PCLU);
    
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
        std::cout<<"before solve"<<std::endl;
        _nl->solve();
        std::cout<<"after solve"<<std::endl;
    }
    
    if (_solve)
    _nl->update();
    
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
    
    std::cout<<"CIAO sono in Residuo"<<std::endl;
    
    NonlinearImplicitSystem * a;
    FEProblemBase::computeResidualSys(a[0],soln,residual);
    
    PetscVector<Number> ciao(_pp_comm);
    
    if (_use_afc && _is_stab_matrix_assembled)
    {
        //_stab_matrix.vector_mult(ciao,soln);
        _stab_matrix.vector_mult_add(residual,soln);
    }
    
    std::cout<<"ho finito Residuo"<<std::endl;
    
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
    
    //m=_stab_matrix.m();
    //n=_stab_matrix.n();
    
    //m_l=_stab_matrix.local_m();
    //n_l=_stab_matrix.local_n();
    
    //rb=_stab_matrix.row_start();
    //re=_stab_matrix.row_stop();
    
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
                
                if (maxEntry>0.0)
                {
                    Real value=-1.0*maxEntry;
                    _stab_matrix.set(row,col,value);
                    _stab_matrix.add(row,row,maxEntry);
                }
            }
        }
        MatRestoreRow(jac_petsc,row,&ncols,&cols,&val);
        MatRestoreRow(jac_tr_petsc,row,&ncols_tr,&cols_tr,&val_tr);
    }
    
    std::cout<<"prima della chiusura\n";
    _stab_matrix.close();
    
    //_stab_matrix.print_matlab("ciao.m");
    
    _is_stab_matrix_assembled=true;
    
}


