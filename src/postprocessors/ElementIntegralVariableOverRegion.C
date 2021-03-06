//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralVariableOverRegion.h"

registerMooseObject("parrotApp", ElementIntegralVariableOverRegion);

template <>
InputParameters
validParams<ElementIntegralVariableOverRegion>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addCoupledVar("variable", "The name of the variable that this object operates on");
  params.addRequiredParam<Real>("region","region");
  return params;
}

ElementIntegralVariableOverRegion::ElementIntegralVariableOverRegion(const InputParameters & parameters) :
ElementIntegralPostprocessor(parameters),
_u( parameters.isParamValid("variable") ? coupledValue("variable") : _zero),
_hasVariable(parameters.isParamValid("variable")),
_regionID(getMaterialProperty<int>("RegionID")),
_region(getParam<Real>("region"))
{}

Real ElementIntegralVariableOverRegion::computeQpIntegral()
{
    if (_region==_regionID[_qp])
    {
        if (_hasVariable)
        {
            return _u[_qp];
        }
        else
        {
            return 1.0;
        }
    }
    else
        return 0.0;
}
