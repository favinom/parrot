#include "FlowAndTransport.h"

// MOOSE includes
//#include "Assembly.h"
//#include "MooseMesh.h"

registerMooseObject("parrotApp", FlowAndTransport);

template <>
InputParameters
validParams<FlowAndTransport>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("phi","phi");
  params.addRequiredParam<Real>("k","k");
  params.addCoupledVar("pressure",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
  return params;
}

FlowAndTransport::FlowAndTransport(const InputParameters &parameters) :
Material(parameters),
_condInput(getParam<Real>("k")),
_poroInput(getParam<Real>("phi")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_poro(declareProperty<Real>("Porosity")),
_K(declareProperty<RealTensorValue>("PermeabilityTensor")),
_U(declareProperty<RealVectorValue>("VelocityVector"))
{}

void
FlowAndTransport::computeQpProperties()
{
    _K[_qp]=RealTensorValue(1.0,0.0,0.0,
                        0.0,1.0,0.0,
                        0.0,0.0,1.0);

    _poro[_qp]=_poroInput;

    _K[_qp]=_condInput*_K[_qp];
    _U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];

}
