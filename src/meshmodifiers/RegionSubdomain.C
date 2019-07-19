//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegionSubdomain.h"
#include "Conversion.h"
#include "MooseMesh.h"


registerMooseObject("MooseApp", RegionSubdomain);

template <>
InputParameters
validParams<RegionSubdomain>()
{
    
    InputParameters params = validParams<MeshModifier>();
    params.addClassDescription("Changes the subdomain ID of elements either (XOR) inside or outside "
                               "the specified box to the specified ID.");
    params.addRequiredParam<SubdomainID>("block_id",
                                         "Subdomain id to set for inside/outside the bounding box");
    params.addParam<SubdomainName>(
                                   "block_name", "Subdomain name to set for inside/outside the bounding box (optional)");
    
    params.addRequiredParam<UserObjectName>("regionUserObject",
                                            "regionUserObject");
    
    
    return params;
}

RegionSubdomain::RegionSubdomain(const InputParameters & parameters) :
MeshModifier(parameters),
_block_id(parameters.get<SubdomainID>("block_id")),
_localFEProblem(getParam<FEProblem *>("_fe_problem")),
_localUserObject(_localFEProblem[0].getUserObjectBase("regionUserObject")),
_regionUserObject(dynamic_cast<RegionUserObject const &>(_localUserObject))
{}

void RegionSubdomain::modify()
{
    // Check that we have access to the mesh
    if (!_mesh_ptr)
        mooseError("_mesh_ptr must be initialized before calling RegionSubdomain::modify()");
    
    // Loop over the elements
    for (const auto & elem : _mesh_ptr->getMesh().active_element_ptr_range())
    {
        Point punto=elem->centroid();
        bool contains = _regionUserObject.isInside(punto);
        if (contains )
        {
            elem->subdomain_id() = _block_id;
        }

        // Assign block name, if provided
        if (isParamValid("block_name"))
            _mesh_ptr->getMesh().subdomain_name(_block_id) = getParam<SubdomainName>("block_name");
    }
}

