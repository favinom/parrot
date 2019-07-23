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

registerMooseObject("parrotApp", ParrotProblem);

template <>
InputParameters
validParams<ParrotProblem>()
{
  InputParameters params = validParams<FEProblem>();
  return params;
}

ParrotProblem::ParrotProblem(const InputParameters & parameters) :
FEProblem(parameters)
{
    _const_jacobian=true;
    PCCreate(PETSC_COMM_WORLD, &_problem_PC);
    PCSetType(_problem_PC,PCLU);
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
}

void
ParrotProblem::stabilize_A_matrix(SparseMatrix<Number> & jacobian)
{
    
    
    _console << "Stabilize A matrix:: begin  "  << std::endl;
    
    SparseMatrix<Number> A_0_t;
    
    jacobian.get_transpose (A_0_t);
    
    S_matrix=A_0_t;
    
    S_matrix.zero();
    
    
    
    
    {
        utopia::each_read(A_0, [&](const utopia::SizeType i, const utopia::SizeType j, double value){
            if(i!=j)
            {
                double value_1 = 1.0 * value;
                
                //utopia::disp(value_1);
                
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
    
    // S_matrix += transpose(S_matrix);
    
    // S_matrix *=0.5;
    
    S_matrix+=S_diag;
    
    _console << "Stabilize A matrix:: end  "  << std::endl;
    
    
}
