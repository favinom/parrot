//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralVolumePostprocessor.h"
#include "FractureUserObject.h"

registerMooseObject("parrotApp", ElementIntegralVolumePostprocessor);

template <>
InputParameters
validParams<ElementIntegralVolumePostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  params.addRequiredParam<int>("fractureRegionId","fractureRegionId");
  return params;
}

ElementIntegralVolumePostprocessor::ElementIntegralVolumePostprocessor(const InputParameters & parameters) :
ElementIntegralPostprocessor(parameters),
_meshModifierName(getParam<std::string>("fractureMeshModifier")), 
_regionId(getParam<int>("fractureRegionId"))
{}

Real
ElementIntegralVolumePostprocessor::computeQpIntegral()
{
  MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
  FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );

  Real bound = 0.0;
  bool check = _fractureUserObject.isInsideRegion(_q_point[_qp], _regionId, bound);
  if(check)
    return 1.0;
  else
    return 0.0;
}
