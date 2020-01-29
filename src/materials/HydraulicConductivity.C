#include "HydraulicConductivity.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseMesh.h"

registerMooseObject("parrotApp", HydraulicConductivity);

template <>
InputParameters
validParams<HydraulicConductivity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("conductivity","conductivity");
  params.addCoupledVar("pressure",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
  params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  return params;
}

HydraulicConductivity::HydraulicConductivity(const InputParameters &parameters) :
Material(parameters),
_cond(getParam<Real>("conductivity")),
_K(declareProperty<RealTensorValue>("PermeabilityTensor")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_U(declareProperty<RealVectorValue>("VelocityVector"))
//_hasMeshModifier(parameters.isParamValid("fractureMeshModifier"))
{
 
  // if (_hasMeshModifier)
  //   {
  //       _meshModifierName=getParam<std::string>("fractureMeshModifier");
  //   }


}

void
HydraulicConductivity::computeQpProperties()
{

    // if (_hasMeshModifier)
    // {
    //   MeshModifier const & _myMeshModifier = _app.getMeshModifier(_meshModifierName.c_str());
    //   RegionUserObject const & regionUserObject=dynamic_cast<RegionUserObject const &>(_myMeshModifier);
    // }
  

  _K[_qp]=RealTensorValue(1.0,0.0,0.0,
                          0.0,1.0,0.0,
                          0.0,0.0,1.0);

  _K[_qp]=_cond*_K[_qp];

    
  _U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];


}
