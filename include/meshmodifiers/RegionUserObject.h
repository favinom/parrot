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
#include "MeshModifier.h"
#include "MooseMesh.h"

// Forward Declarations
class RegionUserObject;

template <>
InputParameters validParams<RegionUserObject>();

/**
 * A user object that runs over all the nodes and does an aggregation
 * step to compute a single value.
 */
class RegionUserObject : public MeshModifier
{
public:
    RegionUserObject(const InputParameters & parameters);
    
    void modify() override;
    
    bool isInside(RealVectorValue const & point, Real bound=0.0) const ;

    std::vector<int> whichIsInside(RealVectorValue const & point,Real bound=0.0) const ;
    
    virtual bool isInsideRegion(RealVectorValue const & point, int region, Real & bound) const = 0;

    
protected:
    
    int _dim;
    
    int _fn;
};
