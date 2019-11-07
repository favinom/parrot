//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Marker.h"

#include "RegionUserObject.h"

class RegionMarker;

template <>
InputParameters validParams<RegionMarker>();

class RegionMarker : public Marker
{
public:
  RegionMarker(const InputParameters & parameters);

protected:
  virtual Marker::MarkerValue computeElementMarker() override;

  std::string _meshModifierName;

  MeshModifier const & _myMeshModifier;
  RegionUserObject const & regionUserObject;

};

