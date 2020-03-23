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
#include <petsc/private/kspimpl.h>
#include "ksp_parrot_impl.h"
#include "StoreOperators.h"

#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

class ParrotProblem;

template <>
InputParameters validParams<ParrotProblem>();

class ParrotProblem : public FEProblem
{
public:
    ParrotProblem(const InputParameters & parameters);
    
    virtual void initialSetup();
    virtual void timestepSetup();
    virtual void solve();
    
    void computeStabilizationMatrix(SparseMatrix<Number> & jacobian);
    
    
    void computeJacobianSys(NonlinearImplicitSystem & /*sys*/,
                            const NumericVector<Number> & soln,
                            SparseMatrix<Number> & jacobian);
    
    void computeResidualSys(NonlinearImplicitSystem & /*sys*/,
                                      const NumericVector<Number> & soln,
                            NumericVector<Number> & residual);

    virtual void update_sol();

    
    KSP_PARROT * _ksp_ptr;
    
    PC _problem_PC;
    int _factorized;
        
    Parallel::Communicator const & _pp_comm;
    PetscMatrix<Number> _stab_matrix;
    
    bool _use_afc;
    bool _is_stab_matrix_assembled;

    std::vector<int> zero_rows;
    StoreOperators * _storeOperatorsUO;
    bool _hasStoreOperatorsUO;

    UserObjectName * userObjectName;

    std::shared_ptr<PetscMatrix<Number>> _poro_lumped;

    std::shared_ptr<PetscVector<Number>> _dirichlet_bc;

    std::shared_ptr<PetscVector<Number>> _value_dirichlet_bc;

    std::shared_ptr<PetscVector<Number>> _sol_vec;
};

