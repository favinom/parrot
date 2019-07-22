//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegionUserObject.h"

template <>
InputParameters
validParams<RegionUserObject>()
{
  InputParameters params = validParams<MeshModifier>();
  return params;
}

RegionUserObject::RegionUserObject(const InputParameters & parameters) :
MeshModifier(parameters)
{}

void RegionUserObject::modify()
{
    if (_dim != _mesh_ptr->dimension() )
    {
        mooseError("In RegionUserObject, the understood _dim and dimension of the provided mesh are different. Most likely, the errore will be in your script.");
    }
};

bool RegionUserObject::isInside(RealVectorValue const & point) const
{
    for (int i=0; i<_fn; ++i)
    {
        if (isInsideRegion(point,i))
            return true;
        
    }
    return false;
};


std::vector<int> RegionUserObject::whichIsInside(RealVectorValue const & point) const
{
    std::vector<int> _whichFrac;
    _whichFrac.clear();
    
    for (int i=0; i<_fn; ++i)
        
    {
        if (isInsideRegion(point,i))
            _whichFrac.push_back(i);
        
    }
    
    return _whichFrac;
};
