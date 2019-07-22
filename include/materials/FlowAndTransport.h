#ifndef HYDRAULICCONDUCTIVITY_H
#define HYDRAULICCONDUCTIVITY_H

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

protected:

    Real const _condInput;
    Real const _poroInput;
    const VariableGradient &_gradP;
  
    MaterialProperty<Real> &_poro;
    MaterialProperty<RealTensorValue> &_K;
    MaterialProperty<RealVectorValue> &_U;

};

#endif
