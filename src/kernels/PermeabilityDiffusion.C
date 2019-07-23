//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PermeabilityDiffusion.h"

registerMooseObject("parrotApp", PermeabilityDiffusion);

template <>
InputParameters
validParams<PermeabilityDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

PermeabilityDiffusion::PermeabilityDiffusion(const InputParameters & parameters) :
Kernel(parameters),
_K(getMaterialProperty<RealTensorValue>("PermeabilityTensor"))
{}

Real
PermeabilityDiffusion::computeQpResidual()
{
  return _grad_u[_qp] * (_K[_qp] *_grad_test[_i][_qp]);
}

Real
PermeabilityDiffusion::computeQpJacobian()
{
  return _grad_phi[_j][_qp] * (_K[_qp] * _grad_test[_i][_qp]);
}
