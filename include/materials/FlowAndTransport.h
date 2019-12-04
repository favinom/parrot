#pragma once

#include "Material.h"

class FlowAndTransport;

template <>
InputParameters validParams<FlowAndTransport>();

class FlowAndTransport : public Material
{
public:
  FlowAndTransport(const InputParameters &parameters);

protected:
  virtual void computeQpProperties();
  //virtual void computeProperties();

protected:

    Real const _condInput;
    Real const _poroInput;
    Real _dim;
    RealTensorValue _id;
    RealVectorValue _u_elem;

    bool const _isPressureValid;
    bool const _conservativeScheme;

    const VariableGradient &_gradP;
  
    MaterialProperty<Real> &_poro;
    MaterialProperty<RealTensorValue> &_K;
    MaterialProperty<RealVectorValue> &_U;

};
