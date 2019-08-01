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

//Forward Declarations
class RegionMaterial;

template<>
InputParameters validParams<RegionMaterial>();

/**
 * Example material class that defines a few properties.
 */
class RegionMaterial : public Material
{
public:
  RegionMaterial(const InputParameters & parameters);
    

protected:
  virtual void computeQpProperties();


    std::vector<int> _whichFrac;
    
    MaterialProperty<int> &  _regionID;
    MaterialProperty<Real> & _regionIDReal;
    MaterialProperty<Real> & _howManyFractures;
    std::string _meshModifierName;
    
    bool _countFractures;
    
    int _elem_counter;
    
};
