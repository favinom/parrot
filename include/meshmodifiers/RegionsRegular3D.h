//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// C++ includes
#include <string>

// MOOSE includes
#include "MinMaxRegion.h"

// Forward Declarations
class RegionsRegular3D;

template <>
InputParameters validParams<RegionsRegular3D>();

/**
 * A user object that runs over all the nodes and does an aggregation
 * step to compute a single value.
 */
class RegionsRegular3D : public MinMaxRegion
{
public:
    RegionsRegular3D(const InputParameters & parameters);
};
