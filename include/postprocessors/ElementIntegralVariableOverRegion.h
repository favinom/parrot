//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"

class ElementIntegralVariableOverRegion;

template <>
InputParameters validParams<ElementIntegralVariableOverRegion>();

class ElementIntegralVariableOverRegion : public ElementIntegralPostprocessor
{
public:
  ElementIntegralVariableOverRegion(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const VariableValue & _u;
  
    bool _hasVariable;
    
  const MaterialProperty<int> & _regionID;
  int _region;
    
};
