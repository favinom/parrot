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

    
    KSP_PARROT * _ksp_ptr;
    
    PC _problem_PC;
    int _factorized;
    
    
    
    Parallel::Communicator const & _pp_comm;
    PetscMatrix<Number> _stab_matrix;
    
    bool _use_afc;
    bool _change_sol;
    std::string _dc_var;
    bool _is_stab_matrix_assembled;
    std::shared_ptr<PetscMatrix<Number>> jmat, smat;

    std::vector<int> _dc_boundary_id;

    std::vector<std::vector<int> > _dc_variables_id;

    std::vector<int> zero_rows;
    // int stabilize_coeffiecient(NumericVector<Number> & vec_solution,
    //                            NumericVector<Number> & ghosted_solution);

    //void outputStep(ExecFlagType type);

    void stabilize_coeffiecient();

    //const StoreOperators & _operator_storage;

    std::string _userobject_name = "operator_userobject_problem";

    bool shouldUpdateSolution();

    bool updateSolution(/*NumericVector<Number> & vec_solution,
                               NumericVector<Number> & ghosted_solution*/);

    void find_boundary(std::vector<int> &zero_rows, 
                       std::vector<int> &_dc_boundary_id);

    void determine_dc_bnd_var_id(const std::vector<std::string> & BC_var);

    std::vector<std::string> split_string(const std::string & s, char delim);
};

