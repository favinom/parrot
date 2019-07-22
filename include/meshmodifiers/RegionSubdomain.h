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
#include "MooseEnum.h"
#include "MeshModifier.h"
#include "FEProblem.h"

#include "RegionUserObject.h"

// Forward declerations
class RegionSubdomain;

template <>
InputParameters validParams<RegionSubdomain>();

namespace libMesh
{
    class BoundingBox;
}

/**
 * MeshModifier for defining a Subdomain inside or outside of a bounding box
 */
class RegionSubdomain : public MeshModifier
{
public:
    /**
     * Class constructor
     * @param parameters The input parameters
     */
    RegionSubdomain(const InputParameters & parameters);
    
    virtual void modify() override;
    
private:
    
    std::string _meshModifierName;
    
    /// Block ID to assign to the region
    SubdomainID _block_id;
    
};
