//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegionMarker.h"

registerMooseObject("parrotApp", RegionMarker);

template <>
InputParameters
validParams<RegionMarker>()
{
    InputParameters params = validParams<Marker>();
    params.addRequiredParam<std::string>("regionMeshModifier","regionMeshModifier");
    //params.addRequiredParam<std::vector<int>>("blockNum", "which block do you mark?");
    return params;
}

RegionMarker::RegionMarker(const InputParameters & parameters) :
Marker(parameters),
_meshModifierName(getParam<std::string>("regionMeshModifier")),
_myMeshModifier( _app.getMeshModifier(_meshModifierName.c_str()) ),
regionUserObject( dynamic_cast<RegionUserObject const &>(_myMeshModifier) )
{
}

Marker::MarkerValue
RegionMarker::computeElementMarker()
{
	Point center=_current_elem[0].centroid();
	Real hmax=_current_elem[0].hmax();
	hmax=hmax/2.0*std::sqrt(2.0);

	if (regionUserObject.isInside(center,hmax))
		return REFINE;
	else
		return DO_NOTHING;
}
