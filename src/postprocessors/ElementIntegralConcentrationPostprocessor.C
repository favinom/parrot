//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralConcentrationPostprocessor.h"
#include "FractureUserObject.h"

registerMooseObject("parrotApp", ElementIntegralConcentrationPostprocessor);

template <>
InputParameters
validParams<ElementIntegralConcentrationPostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  params.addRequiredParam<int>("fractureRegionId","fractureRegionId");
  return params;
}

ElementIntegralConcentrationPostprocessor::ElementIntegralConcentrationPostprocessor(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
   _meshModifierName(getParam<std::string>("fractureMeshModifier")), 
   _regionId(getParam<int>("fractureRegionId")),
   _u(coupledValue("variable"))
{
  addMooseVariableDependency(mooseVariable());
}

Real
ElementIntegralConcentrationPostprocessor::computeQpIntegral()
{
  // std::string _meshModifierName;

  // params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");

  // _meshModifierName(getParam<std::string>("fractureMeshModifier")),

   MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );

   FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
  
   Real bound = 0.0;
   bool check = _fractureUserObject.isInsideRegion(_q_point[_qp], _regionId, bound);

   if(check==true)
       return _u[_qp];
   else
       return 0.0;
}

