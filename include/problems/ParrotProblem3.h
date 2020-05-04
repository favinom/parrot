//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FEProblem.h"
#include "NonlinearSystem.h"

#include "AntidiffusiveFluxes.h"
#include "StoreOperators.h"

#include "ParrotSolver.h"

class ParrotProblem3;

template <>
InputParameters validParams<ParrotProblem3>();

class ParrotProblem3 : public FEProblem
{
public:
    ParrotProblem3(const InputParameters & parameters);
    
    virtual void initialSetup();
    virtual void solve();
    
    void computeStabilizationMatrix(SparseMatrix<Number> & jacobian);
    
    
    void computeJacobianSys(NonlinearImplicitSystem & /*sys*/,
                            const NumericVector<Number> & soln,
                            SparseMatrix<Number> & jacobian);
    
    void computeResidualSys(NonlinearImplicitSystem & /*sys*/,
                            const NumericVector<Number> & soln,
                            NumericVector<Number> & residual);

    // virtual void update_sol();
        
    Parallel::Communicator const & _pp_comm;
    EquationSystems & _equationSystems;
    TransientNonlinearImplicitSystem & _nl_libMesh;
    SparseMatrix <Number> & _mat_SM;
    NumericVector<Number> & _rhs_NV;
    NumericVector<Number> & _sol_NV;
    

    PetscMatrix<Number> _stab_matrix;

    PetscVector<Number> _res_m;
    bool _use_afc;
    bool _is_stab_matrix_assembled;
    bool _is_jac_matrix_assembled;

    std::vector<int> zero_rows;
    StoreOperators * _storeOperatorsUO;
    AntidiffusiveFluxes *_ComputeAF;
    bool _hasStoreOperatorsUO;
    bool _ComputeAntidiffusiveFluxes;

    UserObjectName * userObjectName;
    UserObjectName * userObjectNameFluxes;

    std::shared_ptr<PetscMatrix<Number>> _poro_lumped;

    std::shared_ptr<PetscMatrix<Number>>   _JMatrix;

    std::shared_ptr<PetscVector<Number>> _dirichlet_bc;

    std::shared_ptr<PetscVector<Number>> _value_dirichlet_bc;

    //std::shared_ptr<PetscVector<Number>> _sol_vec;

    ParrotSolver * parrotSolver;
    int _solverType;

};

