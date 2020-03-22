//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IntegralSolutionOverRegionFast.h"

#include "NonlinearSystem.h"

#include "AssembleVolumeVectors.h"

registerMooseObject("parrotApp", IntegralSolutionOverRegionFast);

template <>
InputParameters
validParams<IntegralSolutionOverRegionFast>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<UserObjectName>("VolumeUserObject","The userobject that stores our operators");
  params.addRequiredParam<int>("region","region");
  params.addRequiredParam<bool>("doDomainSize","doDomainSize");
  return params;
}

IntegralSolutionOverRegionFast::IntegralSolutionOverRegionFast(const InputParameters & parameters) :
GeneralPostprocessor(parameters),
_userObjectName(getParam<UserObjectName>("VolumeUserObject")),
_region(getParam<int>("region")),
_doDomainSize(getParam<bool>("doDomainSize"))
{}

Real IntegralSolutionOverRegionFast::getValue()
{
  AssembleVolumeVectors & assembleVolumeVectors=(_fe_problem.getUserObjectTempl<AssembleVolumeVectors>(_userObjectName));

  std::vector< PetscVector<Number> * > const & volumes = assembleVolumeVectors.getVolumeVectors();

  if (_doDomainSize)
  {
    return volumes.at(_region)[0].sum();
  }
  else
  {
    NonlinearSystem & nonlinearSystem=_fe_problem.getNonlinearSystem();
    NumericVector<Number> const * const & solution=nonlinearSystem.currentSolution();

    return solution[0].dot(volumes.at(_region)[0]);
  }
}
