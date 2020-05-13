//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// C++ includes
#include <string>

#include "libmesh/elem.h"

// MOOSE includes
#include "RegionUserObject.h"

// Forward Declarations
class FractureUserObject;

template <>
InputParameters validParams<FractureUserObject>();

/**
 * A user object that runs over all the nodes and does an aggregation
 * step to compute a single value.
 */
class FractureUserObject : public RegionUserObject
{
public:
    FractureUserObject(const InputParameters & parameters);
    
    virtual bool isInsideRegion(RealVectorValue const & point, int const i, Real & bound) const;


    // THESE HAVE BEEN ADDED TO REFINE THE BOUNDARIES
    bool isOnBoundary(Elem & elem) const ;
    bool isOnBoundaryOfRegion(Elem & elem, int region) const;

        
protected:
    
    std::string _fx_string;
    std::string _fy_string;
    std::string _fz_string;
    std::string _fa1_string;
    std::string _fa2_string;
    std::string _fa3_string;
    std::string _fd1_string;
    std::string _fd2_string;
    std::string _fd3_string;
    
    RealVectorValue * _center;
    RealVectorValue * _rotation;
    RealVectorValue * _dimension;
    RealVectorValue ** _n;
    RealVectorValue *_d;

    std::vector< std::vector<RealVectorValue> > _vertex;
    
    void ComputeNormalsFromAngles(RealVectorValue const & angles,
                                  RealVectorValue & n1,
                                  RealVectorValue & n2,
                                  RealVectorValue & n3);
    
    bool isInsideRegion2D(RealVectorValue const & point, int const i, Real & bound) const;
    bool isInsideRegion3D(RealVectorValue const & point, int const i, Real & bound) const;


    // THESE HAVE BEEN ADDED TO REFINE THE BOUNDARIES
    bool isOnBoundaryOfRegion2D(Elem & elem, int region) const;
    bool isOnBoundaryOfRegion3D(Elem & elem, int region) const;

    
};
