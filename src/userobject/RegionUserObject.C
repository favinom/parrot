//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegionUserObject.h"
#include "SubProblem.h"

template <>
InputParameters
validParams<RegionUserObject>()
{
  InputParameters params = validParams<UserObject>();
  return params;
}

RegionUserObject::RegionUserObject(const InputParameters & parameters) :
UserObject(parameters),
_mesh(_subproblem.mesh())
{
    _dim=_mesh.dimension();
}
