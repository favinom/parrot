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

// Forward Declarations
class ElementIntegralVolumePostprocessor;

template <>
InputParameters validParams<ElementIntegralVolumePostprocessor>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class ElementIntegralVolumePostprocessor : public ElementIntegralPostprocessor
{
public:
  ElementIntegralVolumePostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;
 
   std::string _meshModifierName;
   int  _regionId;
};


