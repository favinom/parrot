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
    
    void computeJacobianSys(NonlinearImplicitSystem & /*sys*/,
                            const NumericVector<Number> & soln,
                            SparseMatrix<Number> & jacobian);

    
    KSP_PARROT * _ksp_ptr;
    
    PC _problem_PC;
    int _factorized;
    
    
    SparseMatrix<Number> S_matrix;
    
};

