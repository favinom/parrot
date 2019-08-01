//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegionsRegular3D.h"

registerMooseObject("parrotApp", RegionsRegular3D);

template <>
InputParameters
validParams<RegionsRegular3D>()
{
    InputParameters params = validParams<MinMaxRegion>();
    return params;
}

RegionsRegular3D::RegionsRegular3D(const InputParameters & parameters) :
MinMaxRegion(parameters)
{
    
    _fn = 7+7+8;
    _dim = 3;
    
    // reg 0
    _regionMin.push_back(RealVectorValue(0.0,0.0,0.0) );
    _regionMax.push_back(RealVectorValue(0.5,0.5,0.5) );
    // reg 1
    _regionMin.push_back(RealVectorValue(0.5,0.0,0.0 ) );
    _regionMax.push_back(RealVectorValue(1.0,0.5,0.5 ) );
    // region 2
    _regionMin.push_back(RealVectorValue(0.0,0.5,0.0) );
    _regionMax.push_back(RealVectorValue(0.5,1.0,0.5) );
    // reg 3
    _regionMin.push_back(RealVectorValue(0.5,0.5,0.0) );
    _regionMax.push_back(RealVectorValue(1.0,1.0,0.5) );
    // reg 4
    _regionMin.push_back(RealVectorValue(0.0,0.0,0.5) );
    _regionMax.push_back(RealVectorValue(0.5,0.5,1.0) );
    // reg5
    _regionMin.push_back(RealVectorValue(0.5,0.0,0.5) );
    _regionMax.push_back(RealVectorValue(1.0,0.5,1.0) );
    // reg6
    _regionMin.push_back(RealVectorValue(0.0,0.5,0.5) );
    _regionMax.push_back(RealVectorValue(0.5,1.0,1.0) );
    //reg7
    _regionMin.push_back(RealVectorValue(0.75,0.75,0.75) );
    _regionMax.push_back(RealVectorValue(1.0,1.0,1.0) );
    //reg8
    _regionMin.push_back(RealVectorValue(0.75,0.5  ,0.75) );
    _regionMax.push_back(RealVectorValue(1.0 ,0.75 ,1.0) );
    // reg9
    _regionMin.push_back(RealVectorValue(0.5,0.75,0.75) );
    _regionMax.push_back(RealVectorValue(0.75,1.0,1.0) );
    // reg 10
    _regionMin.push_back(RealVectorValue(0.5,0.5,0.75) );
    _regionMax.push_back(RealVectorValue(0.75,0.75,1.0) );
    //reg 11
    _regionMin.push_back(RealVectorValue(0.75,0.75,0.5) );
    _regionMax.push_back(RealVectorValue(1.0,1.0,0.75) );
    //reg 12
    _regionMin.push_back(RealVectorValue(0.75,0.5,0.5) );
    _regionMax.push_back(RealVectorValue(1.0,0.75,0.75) );
    //reg 13
    _regionMin.push_back(RealVectorValue(0.5,0.75,0.5) );
    _regionMax.push_back(RealVectorValue(0.75,1.0,0.75) );
    // reg 14
    _regionMin.push_back(RealVectorValue(0.5,0.5,0.5) );
    _regionMax.push_back(RealVectorValue(0.625,0.625,0.625) );
    // reg 15
    _regionMin.push_back(RealVectorValue(0.625,0.5,0.5) );
    _regionMax.push_back(RealVectorValue(0.75,0.625,0.625) );
    // reg16
    _regionMin.push_back(RealVectorValue(0.5,0.625,0.5) );
    _regionMax.push_back(RealVectorValue(0.625,0.75,0.625) );
    // reg 17
    _regionMin.push_back(RealVectorValue(0.625,0.625,0.5) );
    _regionMax.push_back(RealVectorValue(0.75,0.75,0.625) );
    // reg 18
    _regionMin.push_back(RealVectorValue(0.5,0.5,0.625) );
    _regionMax.push_back(RealVectorValue(0.625,0.625,0.75) );
    // reg 19
    _regionMin.push_back(RealVectorValue(0.625,0.5,0.625) );
    _regionMax.push_back(RealVectorValue(0.75,0.625,0.75) );
    // reg 20
    _regionMin.push_back(RealVectorValue(0.5,0.625,0.625) );
    _regionMax.push_back(RealVectorValue(0.625,0.75,0.75) );
    // reg 21
    _regionMin.push_back(RealVectorValue(0.625,0.625,0.625) );
    _regionMax.push_back(RealVectorValue(0.75,0.75,0.75) );

}
