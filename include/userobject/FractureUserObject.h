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
    
    virtual bool isInside(RealVectorValue const & point);

    virtual bool isInsideRegion(RealVectorValue const & point, int const i);

    
    virtual std::vector<int> whichIsInside(RealVectorValue const & point);
    
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
    
    void ComputeNormalsFromAngles(RealVectorValue const & angles,
                                                           RealVectorValue & n1,
                                                           RealVectorValue & n2,
                                                           RealVectorValue & n3);
    
    bool isInsideRegion2D(RealVectorValue const & point, int const i);
    bool isInsideRegion3D(RealVectorValue const & point, int const i);
    
};
