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

// Forward declarations
class MyUO;

template <>
InputParameters validParams<MyUO>();

/**
 * Base class for user-specific data
 */
class MyUO : public GeneralUserObject
{
protected:
    
    FEProblem * _fe_problem;
    EquationSystems & _equationSystems;
    
    
public:
  MyUO(const InputParameters & params);

    virtual void execute() override ;

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize() override ;

  /**
   * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where
   * you want to do MPI communication!
   */
  virtual void finalize() override ;

  /**
   * Must override.
   *
   * @param uo The UserObject to be joined into _this_ object.  Take the data from the uo object and
   * "add" it into the data for this object.
   */
    virtual void threadJoin(const UserObject & uo) override;

};

