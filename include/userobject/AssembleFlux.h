//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
//#include "UserObject.h"

#include "GeneralUserObject.h"
#include "MeshModifier.h"
#include "FractureUserObject.h"

#include "libmesh/quadrature.h"

// Forward declarations
class AssembleFlux;

template <>
InputParameters validParams<AssembleFlux>();

/**
 * Base class for user-specific data

 */

class AssembleFlux : public GeneralUserObject
{

    
    
public:
  AssembleFlux(const InputParameters & params);

  virtual void execute() override {};

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize() override;

  /**
   * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where
   * you want to do MPI communication!
   */
  virtual void finalize() override {};

  /**
   * Must override.
   *
   * @param uo The UserObject to be joined into _this_ object.  Take the data from the uo object and
   * "add" it into the data for this object.
   */
  virtual void threadJoin(const UserObject &) override {};

  protected:

//    Parallel::Communicator const & _pp_comm;
    
//    PetscMatrix<Number> _stiffness_matrix_1;
//    PetscMatrix<Number> _stiffness_matrix_2;
//    PetscMatrix<Number> _stiffness_matrix_t;
    
    AuxVariableName _sol_var_name;

    std::vector<int> _vector_p;
    std::vector<Real> _vector_value;
    std::vector<boundary_id_type> _boundary_D_ids;
    std::vector<boundary_id_type> _boundary_N_ids;
    std::vector<Real> _value_N_bc;
    std::vector<Real> _value_D_bc;

    void ComputeFlux();

    Real ComputeMaterialProprties(const Elem *elem);

    unsigned int _p_var;

    std::string _var_name;
    std::string _sys_name;

    int solve();


    bool _hasMeshModifier;
    std::string _meshModifierName;

    bool const _conservativeScheme;

    QBase const * const & _qrule;


};

