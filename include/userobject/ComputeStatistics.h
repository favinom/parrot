//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

// Forward declarations
class ComputeStatistics;

template <>
InputParameters validParams<ComputeStatistics>();

class ComputeStatistics : public GeneralUserObject
{
public:
 ComputeStatistics(const InputParameters & params);

  virtual void initialize() override {};
  virtual void execute() override ;
  virtual void finalize() override {};
  virtual void threadJoin(const UserObject & ) override {};

  protected:

    void assemble_mass_matrix();

    bool _hasMeshModifier;
    std::string _meshModifierName;

};

