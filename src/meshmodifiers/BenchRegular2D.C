//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BenchRegular2D.h"

registerMooseObject("parrotApp", BenchRegular2D);

template <>
InputParameters
validParams<BenchRegular2D>()
{
    InputParameters params = validParams<MinMaxRegion>();
    return params;
}

BenchRegular2D::BenchRegular2D(const InputParameters & parameters) :
MinMaxRegion(parameters)
{
    
    _fn = 10;
    _dim = 2;
    
    _regionMin.push_back(RealVectorValue(0.0,0.0,0.0) );
    _regionMax.push_back(RealVectorValue(0.5,0.5,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.5,0.0,0.0) );
    _regionMax.push_back(RealVectorValue(1.01,0.5,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.0,0.5,0.0) );
    _regionMax.push_back(RealVectorValue(0.5,1.01,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.750,0.750,0.0) );
    _regionMax.push_back(RealVectorValue(1.01,1.01,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.75,0.50,0.0) );
    _regionMax.push_back(RealVectorValue(1.01,0.75,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.5,0.750,0.0) );
    _regionMax.push_back(RealVectorValue(0.75,1.01,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.50,0.50,0.0) );
    _regionMax.push_back(RealVectorValue(0.625,0.625,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.6250,0.50,0.0) );
    _regionMax.push_back(RealVectorValue(0.75,0.625,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.5,0.6250,0.0) );
    _regionMax.push_back(RealVectorValue(0.625,0.75,0.0) );
    
    _regionMin.push_back(RealVectorValue(0.625,0.6250,0.0) );
    _regionMax.push_back(RealVectorValue(0.75,0.75,0.0) );
    
}
