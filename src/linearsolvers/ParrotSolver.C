#include"ParrotSolver.h"

ParrotSolver::ParrotSolver(int solverType) :
ParrotSolver(solverType,PETSC_COMM_WORLD)
{}

ParrotSolver::ParrotSolver(int solverType, Parallel::Communicator const & pp_comm) :
ParrotSolver(solverType,pp_comm.get())
{}

ParrotSolver::ParrotSolver(int solverType, MPI_Comm const & comm) : _solverType(solverType), _communicator(comm), _constantMatrix(false), _factorized(false)
{
  Construct();
}

PetscErrorCode ParrotSolver::Construct()
{
  _ierr = KSPCreate(_communicator,&_ksp); CHKERRQ(_ierr);
  if (_solverType == 1 || _solverType == 2 )
  {
    _ierr = KSPSetType(_ksp,KSPPREONLY); CHKERRQ(_ierr);
  }
  if (_solverType == 2 )
  {
    _ierr = KSPSetType(_ksp,KSPGMRES); CHKERRQ(_ierr); // KSPRICHARDSON
  }

  _ierr = KSPGetPC(_ksp,&_pc); CHKERRQ(_ierr);

  if (_solverType == 1)
  {
    _ierr = PCSetType(_pc,PCLU); CHKERRQ(_ierr);
    //_ierr = PCFactorSetMatSolverType(_pc,MATSOLVERMUMPS); CHKERRQ(_ierr);
    _ierr = PCFactorSetMatSolverPackage(_pc,MATSOLVERMUMPS); CHKERRQ(_ierr);
    
  }
  if (_solverType == 2)
  {
    _ierr = PCSetType(_pc,PCCHOLESKY); CHKERRQ(_ierr);
    //_ierr = PCFactorSetMatSolverType(_pc,MATSOLVERMUMPS); CHKERRQ(_ierr);
    _ierr = PCFactorSetMatSolverPackage(_pc,MATSOLVERMUMPS); CHKERRQ(_ierr);
  }
  if (_solverType == 3)
  {
    _ierr = PCSetType(_pc,PCHYPRE); CHKERRQ(_ierr);
    _ierr = PCHYPRESetType(_pc, "boomeramg"); CHKERRQ(_ierr);
  }
  PetscFunctionReturn(0);
}

ParrotSolver::~ParrotSolver()
{
  KSPDestroy(&_ksp);
  //PCDestroy(&_pc);
  VecDestroy(&_residual);
  MatDestroy(&_preconditioner);

}

PetscErrorCode ParrotSolver::setMatrixAndVectors(
  SparseMatrix <Number> * & mat_SM,
  NumericVector<Number> * & rhs_NV,
  NumericVector<Number> * & sol_NV)
{
  _sol_NV=sol_NV;
  _rhs_NV=rhs_NV;
  _mat_SM=mat_SM;

  _sol_PV=dynamic_cast<PetscVector<Number> *>(_sol_NV);
  _rhs_PV=dynamic_cast<PetscVector<Number> *>(_rhs_NV);
  _mat_PM=dynamic_cast<PetscMatrix<Number> *>(_mat_SM);

  _ierr = VecDuplicate(_rhs_PV->vec(),&_residual);
  _ierr = MatDuplicate(_mat_PM->mat(),MAT_COPY_VALUES,&_preconditioner);

  //_ierr = PCSetOperators(_pc, _mat_PM->mat(),_preconditioner); CHKERRQ(_ierr);
  _ierr = KSPSetOperators(_ksp, _mat_PM->mat(),_preconditioner); CHKERRQ(_ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ParrotSolver::solve()
{     
  if (_factorized==false || _constantMatrix==false)
  {
    std::cout<<"Factorizing"<<std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();
    //_ierr = PCSetUp(_pc); CHKERRQ(_ierr);
    _ierr = KSPSetUp(_ksp); CHKERRQ(_ierr);
    auto t_stop = std::chrono::high_resolution_clock::now();
    _factorized=true;

    std::cout<<"Factorization time: "<< std::chrono::duration<double, std::milli>(t_stop-t_start).count()<< " ms\n";

  }

  if (_constantMatrix==true)
  {
    //_ierr = PCSetReusePreconditioner(_ksp, PETSC_TRUE); CHKERRQ(_ierr);
    _ierr = KSPSetReusePreconditioner(_ksp, PETSC_TRUE); CHKERRQ(_ierr);
  }

  std::cout<<"Solving"<<std::endl;
  computeResidual();
  auto t_start = std::chrono::high_resolution_clock::now();
  //_ierr = PCApply(_pc,_rhs_PV->vec(),_sol_PV->vec()); CHKERRQ(_ierr);
  _ierr = KSPSolve(_ksp,_rhs_PV->vec(),_sol_PV->vec()); CHKERRQ(_ierr);
  auto t_stop = std::chrono::high_resolution_clock::now();

  std::cout<<"Solving time: "<< std::chrono::duration<double, std::milli>(t_stop-t_start).count()<< " ms\n";
  computeResidual();

  PetscFunctionReturn(0);

}

PetscErrorCode ParrotSolver::computeResidual()
{

  PetscReal norm;
  _ierr =  MatResidual(_mat_PM->mat(),_rhs_PV->vec(),_sol_PV->vec(),_residual); CHKERRQ(_ierr);
  _ierr =  VecNorm(_residual,NORM_2,&norm); CHKERRQ(_ierr);
  std::cout<<"Residual norm: "<<norm<<std::endl;
  PetscFunctionReturn(0);
}

//
      // std::cout<<"start solving?\n";
      // t_start = std::chrono::high_resolution_clock::now();

    
//     if (factorized==0)
//     {
//         _ksp_ptr[0].factorized[0]=1;
    
    
        
//     std::cout<<"start factorizing?\n";
//     auto t_start = std::chrono::high_resolution_clock::now();
//     
//     auto t_end = std::chrono::high_resolution_clock::now();
//     std::cout<<"done factorizing?\n";
//         std::cout<<"fact time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";
//     }
//     PCSetReusePreconditioner(_ksp_ptr[0].local_pc[0], PETSC_TRUE);
//     std::cout<<"start solving?\n";
//     auto t_start = std::chrono::high_resolution_clock::now();
//     PCApply(_ksp_ptr[0].local_pc[0],ksp->vec_rhs,ksp->vec_sol);
//     auto t_end = std::chrono::high_resolution_clock::now();
//     std::cout<<"done solving?\n";
//     std::cout<<"solve time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";

//     _ksp_ptr[0].local_pc=NULL;
    
// //    std::cout<<"start solving?\n";
// //     t_start = std::chrono::high_resolution_clock::now();
// //    PCApply(_ksp_ptr[0].local_pc[0],ksp->vec_rhs,ksp->vec_sol);
// //     t_end = std::chrono::high_resolution_clock::now();
// //    std::cout<<"done solving?\n";
// //    std::cout<<"solve time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";
    
//     Vec r;
//     VecDuplicate(ksp->vec_rhs,&r);
//     MatResidual(Hmat,ksp->vec_rhs,ksp->vec_sol,r);
//     PetscReal norm;
//     VecNorm(r,NORM_2,&norm);
//     std::cout<<"qui "<<norm<<std::endl;
//     PetscPrintf(PETSC_COMM_WORLD,"   %14.12e \n", norm);

    
// //  ierr     = KSP_PCApply(ksp,ksp->vec_rhs,ksp->vec_sol);CHKERRQ(ierr);
// //  ierr     = PCGetSetUpFailedReason(ksp->pc,&pcreason);CHKERRQ(ierr);
// //  if (pcreason) {
// //    ksp->reason = KSP_DIVERGED_PCSETUP_FAILED;
// //  } else {
//     ksp->its    = 1;
//     ksp->reason = KSP_CONVERGED_ITS;
// //  }

