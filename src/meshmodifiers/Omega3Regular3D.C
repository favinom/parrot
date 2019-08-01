//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Omega3Regular3D.h"

registerMooseObject("parrotApp", Omega3Regular3D);

template <>
InputParameters
validParams<Omega3Regular3D>()
{
    InputParameters params = validParams<MinMaxRegion>();
    return params;
}

Omega3Regular3D::Omega3Regular3D(const InputParameters & parameters) :
MinMaxRegion(parameters)
{
    
    _fn = 3;
    _dim = 3;
    
    _regionMin.push_back(RealVectorValue(0.5,0.0,0.0) );
    _regionMax.push_back(RealVectorValue(1.0,0.5,1.0) );
    
    _regionMin.push_back(RealVectorValue(0.75,0.5,0.5) );
    _regionMax.push_back(RealVectorValue(1.01,0.75,1.0) );
    
    _regionMin.push_back(RealVectorValue(0.625,0.5,0.5) );
    _regionMax.push_back(RealVectorValue(0.75,0.625,0.75) );
    
}
