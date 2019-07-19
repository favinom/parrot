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
#include "UserObject.h"

// Forward Declarations
class RegionUserObject;

template <>
InputParameters validParams<RegionUserObject>();

/**
 * A user object that runs over all the nodes and does an aggregation
 * step to compute a single value.
 */
class RegionUserObject : public UserObject
{
public:
  RegionUserObject(const InputParameters & parameters);
    
    virtual bool isInside(RealVectorValue const & point) = 0;

    virtual bool isInsideRegion(RealVectorValue const & point, int region) = 0;
    
    virtual std::vector<int> whichIsInside(RealVectorValue const & point) = 0;
    
protected:
  /// The mesh that is being iterated over
  MooseMesh & _mesh;

    int _dim;
    
    int _fn;
};
