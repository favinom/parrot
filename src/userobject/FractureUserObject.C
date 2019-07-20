//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FractureUserObject.h"
#include "SubProblem.h"

registerMooseObject("parrotApp", FractureUserObject);

template <>
InputParameters
validParams<FractureUserObject>()
{
    InputParameters params = validParams<RegionUserObject>();
    
    params.addRequiredParam<int>("fn", "number of fractures");
    params.addRequiredParam<std::string>("fx_string", "x-coordinates of center of fractures");
    params.addRequiredParam<std::string>("fy_string", "y-coordinates of center of fractures");
    params.addParam<std::string>("fz_string", "z-coordinates of center of fractures");
    params.addRequiredParam<std::string>("fa1_string", "rotation along z axis");
    params.addParam<std::string>("fa2_string", "rotation along y axis");
    params.addParam<std::string>("fa3_string", "rotation along x axis");
    params.addRequiredParam<std::string>("fd1_string", "fracture dimension 1");
    params.addRequiredParam<std::string>("fd2_string", "fracture dimension 2");
    params.addParam<std::string>("fd3_string", "fracture dimension 3");
    
    
    return params;
}

FractureUserObject::FractureUserObject(const InputParameters & parameters) :
RegionUserObject(parameters),
_fx_string(getParam<std::string>("fx_string")),
_fy_string(getParam<std::string>("fy_string")),
_fa1_string(getParam<std::string>("fa1_string")),
_fd1_string(getParam<std::string>("fd1_string")),
_fd2_string(getParam<std::string>("fd2_string"))
{
    MeshModifier const & _myMeshModifier= _app.getMeshModifier("ciao");
    
    _fn=getParam<int>("fn");
    
    if (_dim==2 && parameters.isParamValid("fz_string") )
        mooseError("You provided fz_string in a 2D code");
    
    if (_dim==2 && parameters.isParamValid("fa2_string") )
        mooseError("You provided fa2_string in a 2D code");
    
    if (_dim==2 && parameters.isParamValid("fa3_string") )
        mooseError("You provided fa3_string in a 2D code");
    
    if (_dim==2 && parameters.isParamValid("fd3_string") )
        mooseError("You provided fd3_string in a 2D code");
    
    if (_dim==3 && !parameters.isParamValid("fz_string") )
        mooseError("You did not provided fz_string in a 3D code");
    
    if (_dim==3 && !parameters.isParamValid("fa2_string") )
        mooseError("You did not provided fa2_string in a 3D code");
    
    if (_dim==3 && !parameters.isParamValid("fa3_string") )
        mooseError("You did not provided fa3_string in a 3D code");
    
    if (_dim==3 && !parameters.isParamValid("fd3_string") )
        mooseError("You did not provided fd3_string in a 3D code");
    
    
    if (_dim==3)
    {
        _fz_string=getParam<std::string>("fz_string");
        _fa2_string=getParam<std::string>("fa2_string");
        _fa3_string=getParam<std::string>("fa3_string");
        _fd3_string=getParam<std::string>("fd3_string");
    }
    
    _center   =new RealVectorValue [_fn];
    _rotation =new RealVectorValue [_fn];
    _dimension=new RealVectorValue [_fn];
    _d        =new RealVectorValue [_fn];
    
    _n = new RealVectorValue * [_fn];
    for (int i=0; i<_fn; ++i)
    {
        _n[i]=new RealVectorValue [_dim];
    }
    
    
    std::istringstream  fx_ss( _fx_string);
    std::istringstream  fy_ss( _fy_string);
    std::istringstream fa1_ss(_fa1_string);
    std::istringstream fd1_ss(_fd1_string);
    std::istringstream fd2_ss(_fd2_string);
    
    std::string token;
    
    for (int i=0; i<_fn; ++i)
    {
        std::getline(fx_ss, token, ',');
        _center[i](0)=std::atof(token.c_str());
        std::getline(fy_ss, token, ',');
        _center[i](1)=std::atof(token.c_str());
        
        std::getline(fa1_ss, token, ',');
        _rotation[i](0)=atof(token.c_str())/180.0*pi;
        _rotation[i](1)=0.0;
        _rotation[i](2)=0.0;
        
        std::getline(fd1_ss, token, ',');
        _dimension[i](0)=atof(token.c_str());
        std::getline(fd2_ss, token, ',');
        _dimension[i](1)=atof(token.c_str());
        
        for (int j=0; j<2; ++j)
        {
            _d[i](j)=_n[i][j]*_center[i];
        }
    }
    
    
    if (_dim==3)
    {
        
        
        std::istringstream  fz_ss( _fz_string);
        std::istringstream fa2_ss(_fa2_string);
        std::istringstream fa3_ss(_fa3_string);
        std::istringstream fd3_ss(_fd3_string);
        
        for (int i=0; i<_fn; ++i)
        {
            std::getline(fz_ss, token, ',');
            _center[i](2)=std::atof(token.c_str());
            
            std::getline(fa2_ss, token, ',');
            _rotation[i](1)=atof(token.c_str())/180.0*pi;
            std::getline(fa3_ss, token, ',');
            _rotation[i](2)=atof(token.c_str())/180.0*pi;
            
            std::getline(fd3_ss, token, ',');
            _dimension[i](2)=atof(token.c_str());
            
            
            _d[i](3)=_n[i][3]*_center[i];
        }
        
    }
    
    for (int i=0; i<_fn; ++i)
    {
        ComputeNormalsFromAngles(_rotation[i],_n[i][0],_n[i][1],_n[i][2]);
    }
    
    delete [] _center;
    delete [] _rotation;
    
    std::cout<<"costruito lo UO\n";
    
}

bool FractureUserObject::isInsideRegion(RealVectorValue const & point, int const i) const
{
    bool ret;
    if (_dim==2)
    {
        ret=isInsideRegion2D(point,i);
    }
    if (_dim==3)
    {
        ret=isInsideRegion3D(point,i);
    }
    return ret;
};

bool FractureUserObject::isInside(RealVectorValue const & point) const
{
    for (int i=0; i<_fn; ++i)
    {
        if (isInsideRegion(point,i))
            return true;
        
    }
    return false;
};



std::vector<int> FractureUserObject::whichIsInside(RealVectorValue const & point) const
{
    std::vector<int> _whichFrac;
    _whichFrac.clear();
    
    for (int i=0; i<_fn; ++i)
    
    {
        if (isInsideRegion(point,i))
            _whichFrac.push_back(i);
        
    }
    
    return _whichFrac;
};


void FractureUserObject::ComputeNormalsFromAngles(RealVectorValue const & angles,
                                                  RealVectorValue & n1,
                                                  RealVectorValue & n2,
                                                  RealVectorValue & n3)
{
    RealTensorValue R1;
    RealTensorValue R2;
    RealTensorValue R3;
    
    R1(0,0)=std::cos(angles(0));
    R1(0,1)=-std::sin(angles(0));
    R1(0,2)=0.0;
    R1(1,0)=std::sin(angles(0));
    R1(1,1)=std::cos(angles(0));
    R1(1,2)=0.0;
    R1(2,0)=0.0;
    R1(2,1)=0.0;
    R1(2,2)=1.0;
    
    R2(0,0)=std::cos(angles(1));
    R2(0,1)=0.0;
    R2(0,2)=-std::sin(angles(1));
    R2(1,0)=0.0;
    R2(1,1)=1.0;
    R2(1,2)=0.0;
    R2(2,0)=std::sin(angles(1));
    R2(2,1)=0.0;
    R2(2,2)=std::cos(angles(1));
    
    R3(0,0)=1.0;
    R3(0,1)=0.0;
    R3(0,2)=0.0;
    R3(1,0)=0.0;
    R3(1,1)=std::cos(angles(2));
    R3(1,2)=-std::sin(angles(2));
    R3(2,0)=0.0;
    R3(2,1)=std::sin(angles(2));
    R3(2,2)=std::cos(angles(2));
    
    RealTensorValue R=R1*R2*R3;
    
    for (int i=0; i<3; ++i)
    {
        n1(i)=R(i,0);
        n2(i)=R(i,1);
        n3(i)=R(i,2);
    }
}


bool FractureUserObject::isInsideRegion2D(RealVectorValue const & point, int const i) const
{
    Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
    if (temp1<_dimension[i](0)/2.0)
    {
        Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
        if (temp2<_dimension[i](1)/2.0)
        {
            return true;
            
        }
    }
    return false;
    
}

bool FractureUserObject::isInsideRegion3D(RealVectorValue const & point, int const i) const
{
    Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
    if (temp1<_dimension[i](0)/2.0)
    {
        Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
        if (temp2<_dimension[i](1)/2.0)
        {
            Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
            if (temp3<_dimension[i](2)/2.0)
            {
                return true;
            }
            
        }
    }
    return false;
}
