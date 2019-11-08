//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MinMaxRegion.h"
#include "MooseApp.h"

registerMooseObject("parrotApp", MinMaxRegion);

template <>
InputParameters
validParams<MinMaxRegion>()
{
    InputParameters params = validParams<RegionUserObject>();
    return params;
}

MinMaxRegion::MinMaxRegion(const InputParameters & parameters) :
RegionUserObject(parameters)
{}

bool MinMaxRegion::isInsideRegion(RealVectorValue const & point, int const i, Real & bound) const
{
    bool ret;
    if (_dim==2)
    {
        ret=isInsideRegion2D(point,i,bound);
    }
    if (_dim==3)
    {
        ret=isInsideRegion3D(point,i,bound);
    }
    return ret;
};


bool MinMaxRegion::isInsideRegion2D(RealVectorValue const & point, int const i, Real & bound) const
{
    bool isIn=true;
    for (int dim=0; dim<2; ++dim)
    {
        if ( _regionMin[i](dim) < point(dim) && point(dim) < _regionMax[i](dim) )
        {
            // do nothing
        }
        else
        {
            isIn = false;
        }
    }

    return isIn;
}
    
bool MinMaxRegion::isInsideRegion3D(RealVectorValue const & point, int const i, Real & bound) const
{
    bool isIn=true;
    for (int dim=0; dim<3; ++dim)
    {
        if ( _regionMin[i](dim) < point(dim) && point(dim) < _regionMax[i](dim) )
        {
            // do nothing
        }
        else
        {
            isIn = false;
        }
    }
    
    return isIn;
}
