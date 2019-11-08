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
#include "RegionUserObject.h"

// Forward Declarations
class MinMaxRegion;

template <>
InputParameters validParams<MinMaxRegion>();

/**
 * A user object that runs over all the nodes and does an aggregation
 * step to compute a single value.
 */
class MinMaxRegion : public RegionUserObject
{
public:
    MinMaxRegion(const InputParameters & parameters);
    
    virtual bool isInsideRegion(RealVectorValue const & point, int const i, Real & bound) const;
        
protected:
    
    
    bool isInsideRegion2D(RealVectorValue const & point, int const i, Real & bound) const;
    bool isInsideRegion3D(RealVectorValue const & point, int const i, Real & bound) const;
    
    std::vector<RealVectorValue> _regionMin;
    std::vector<RealVectorValue> _regionMax;
    
};
