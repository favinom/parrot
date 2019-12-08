/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "FractureMaterial.h"

registerMooseObject("parrotApp", FractureMaterial);

template<>
InputParameters validParams<FractureMaterial>()
{
    InputParameters params = validParams<Material>();
    params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");
    params.addRequiredParam<Real>("matrixPorosity","matrixPorosity");
    params.addRequiredParam<Real>("fracturePorosity","fracturePorosity");
    params.addRequiredParam<Real>("matrixPermeability","matrixPermeability");
    params.addRequiredParam<Real>("fracturePermeability","fracturePermeability");
    params.addCoupledVar("pressure",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
    
    return params;
}

FractureMaterial::FractureMaterial(const InputParameters & parameters) :
Material(parameters),
_meshModifierName(getParam<std::string>("fractureMeshModifier")),
//_myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) ),
//_fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) ),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_poro(declareProperty<Real>("Porosity")),
_K(declareProperty<RealTensorValue>("PermeabilityTensor")),
_U(declareProperty<RealVectorValue>("VelocityVector")),
_poroMatrix(getParam<Real>("matrixPorosity")),
_poroFracture(getParam<Real>("fracturePorosity")),
_kappaMatrix(getParam<Real>("matrixPermeability")),
_kappaFracture(getParam<Real>("fracturePermeability"))
{
    //_poroMatrix=1.0;
    //_poroFracture=1.0;
    //_kappaMatrix=1e-8;
    //_kappaFracture=1e-6;
}

void
FractureMaterial::computeQpProperties()
{ 

    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );

    _K[_qp]=RealTensorValue(1.0,0.0,0.0,
                        0.0,1.0,0.0,
                        0.0,0.0,1.0);

    if ( _fractureUserObject.isInside(_q_point[_qp]) )
    {
        _poro[_qp]=_poroFracture;
        _K[_qp]=_kappaFracture*_K[_qp];
    }
    else
    {
        _poro[_qp]=_poroMatrix;
        _K[_qp]=_kappaMatrix*_K[_qp];
    }

    _U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];
    
}
