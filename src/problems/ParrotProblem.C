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
    params.addParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    return params;
}

ParrotProblem::ParrotProblem(const InputParameters & parameters) :
FEProblem(parameters),
_pp_comm(_mesh.comm()),
_stab_matrix(_pp_comm),
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
        
    NonlinearImplicitSystem * a;
    FEProblemBase::computeResidualSys(a[0],soln,residual);
    
    if (_use_afc && _is_stab_matrix_assembled)
    {
        _stab_matrix.vector_mult_add(residual,soln);
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

    //_stab_matrix = const_cast<StoreOperators&>(getUserObjectTempl<StoreOperators>(_userobject_name)).StabMatrix();
    
    // initialize stabilization matrix
    _stab_matrix.init(m,n,m_l,n_l,30);
    
    Mat _stab_matrix_petsc=_stab_matrix.mat();
    MatSetOption(_stab_matrix_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    
    MatCopy(jac_PM.mat(),  _stab_matrix.mat(), DIFFERENT_NONZERO_PATTERN);
    _stab_matrix.zero();

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
        }
        MatRestoreRow(jac_petsc,row,&ncols,&cols,&val);
        MatRestoreRow(jac_tr_petsc,row,&ncols_tr,&cols_tr,&val_tr);
    }
    
     _stab_matrix.close();
        
    _is_stab_matrix_assembled=true;

    std::cout<<"STOP ParrotProblem::computeStabilizationMatrix\n";
}
