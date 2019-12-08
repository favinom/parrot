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

#pragma once

#include "Material.h"
#include "FractureUserObject.h"

//Forward Declarations
class FractureMaterial;

template<>
InputParameters validParams<FractureMaterial>();

/**
 * Example material class that defines a few properties.
 */
class FractureMaterial : public Material
{
public:
  FractureMaterial(const InputParameters & parameters);
    

protected:
  virtual void computeQpProperties();


    //std::vector<int> _whichFrac;
    
    //MaterialProperty<int> &  _regionID;
    //MaterialProperty<Real> & _regionIDReal;
    //MaterialProperty<Real> & _howManyFractures;
    std::string _meshModifierName;

//    MeshModifier const & _myMeshModifier;
//    FractureUserObject const & _fractureUserObject;
    
    const VariableGradient &_gradP;

    MaterialProperty<Real> &_poro;
    MaterialProperty<RealTensorValue> &_K;
    MaterialProperty<RealVectorValue> &_U;


    Real const _poroMatrix;
    Real const _poroFracture;
    Real const _kappaMatrix;
    Real const _kappaFracture;
    //bool _countFractures;
    
    //int _elem_counter;
    
};
