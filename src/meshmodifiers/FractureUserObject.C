//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FractureUserObject.h"
#include "MooseApp.h"

#include "geohelpers.h"

#define myeps 1e-15

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
    // we read the number of fractures
    _fn=getParam<int>("fn");
    // if all non required params are valid we assume dimension is 3, otherwise dimension is 2
    if (parameters.isParamValid("fz_string")  &&
        parameters.isParamValid("fa2_string") &&
        parameters.isParamValid("fa3_string") &&
        parameters.isParamValid("fd3_string")   )
    {
        
        std::cout<<"Assuming dimension is 3\n";
        _dim=3;
    }
    else
    {
        std::cout<<"Assuming dimension is 2\n";
        _dim=2;
    }
    
    // Here we check inconsistencies:
    // if _dim=2 but any non-required parameter is provided something is wrong
    if (_dim==2 &&  parameters.isParamValid("fz_string") )
        mooseError("You provided fz_string in a 2D code");
    
    if (_dim==2 && parameters.isParamValid("fa2_string") )
        mooseError("You provided fa2_string in a 2D code");
    
    if (_dim==2 && parameters.isParamValid("fa3_string") )
        mooseError("You provided fa3_string in a 2D code");
    
    if (_dim==2 && parameters.isParamValid("fd3_string") )
        mooseError("You provided fd3_string in a 2D code");
    
    // if dim is 3, we read the non-required params
    if (_dim==3)
    {
        _fz_string=getParam<std::string>("fz_string");
        _fa2_string=getParam<std::string>("fa2_string");
        _fa3_string=getParam<std::string>("fa3_string");
        _fd3_string=getParam<std::string>("fd3_string");
    }
    
    // At this point, we have read all parameters and set the dimension
    
    // we allocate the arrays to store the parameters
    
    // coordinates of the centers
    _center   =new RealVectorValue [_fn];
    // rotation angles
    _rotation =new RealVectorValue [_fn];
    // dimension
    _dimension=new RealVectorValue [_fn];
    // the _d defining the planes
    _d        =new RealVectorValue [_fn];
    
    // the normals and vertices
    _n = new RealVectorValue * [_fn];
    _vertex.resize(_fn);
    for (int i=0; i<_fn; ++i)
    {
        // it is correct that we allocate with dimension 3, it is necessary when we call compute normals from angles
        _n[i]=new RealVectorValue [3];
        if (_dim==2)
            _vertex.at(i).resize(4);
        else
            if (_dim==3)
                _vertex.at(i).resize(8);
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
        _center[i](2)=0.0;
        
        std::getline(fa1_ss, token, ',');
        _rotation[i](0)=atof(token.c_str())/180.0*pi;
        _rotation[i](1)=0.0;
        _rotation[i](2)=0.0;
        
        std::getline(fd1_ss, token, ',');
        _dimension[i](0)=atof(token.c_str());
        std::getline(fd2_ss, token, ',');
        _dimension[i](1)=atof(token.c_str());
        _dimension[i](2)=0.0;
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
            
        }
        
    }

    for (int i=0; i<_fn; ++i)
    {
        // maybe one day, we write two different functions: 2D and 3D
        ComputeNormalsFromAngles(_rotation[i],_n[i][0],_n[i][1],_n[i][2]);
        
        for (int j=0; j<2; ++j)
        {
            _d[i](j)=_n[i][j]*_center[i];
        }
        
        if (_dim==3)
            _d[i](2)=_n[i][2]*_center[i];
        
    }

    for (int i=0; i<_fn; ++i)
    {
        RealTensorValue R;
        for (int l=0; l<3; ++l)
        {
            R(l,0)=_n[i][0](l);
            R(l,1)=_n[i][1](l);
            R(l,2)=_n[i][2](l);
        }
        int counter=0;

        if (_dim==3)
        {
            for (int kk=-1; kk<=1; ++++kk)
                for (int jj=-1; jj<=1; ++++jj)
                    for (int ii=-1; ii<=1; ++++ii)
                    {
                        RealVectorValue & t=_vertex.at(i).at(counter);
                        t=RealVectorValue( ii*0.5*_dimension[i](0),jj*0.5*_dimension[i](1),kk*0.5*_dimension[i](2) );
                        t=R*t;
                        t=t+_center[i];
                        counter++;
                    }
        }
        else
            if (_dim==2)
            {
                for (int jj=-1; jj<=1; ++++jj)
                    for (int ii=-1; ii<=1; ++++ii)
                    {
                        RealVectorValue & t=_vertex.at(i).at(counter);
                        t=RealVectorValue( ii*0.5*_dimension[i](0),jj*0.5*_dimension[i](1),0.0 );
                        t=R*t;
                        t=t+_center[i];
                        counter++;
                    }
            }
    }

    delete [] _center;
    delete [] _rotation;

}

bool FractureUserObject::isInsideRegion(RealVectorValue const & point, int const i, Real & bound) const
{
    bool ret;
    if (_dim==2)
    {
        ret=isInsideRegion2D(point,i,bound);
    }
    if (_dim==3)
    {
        ret=isInsideRegion3D(point,i,bound);
    }
    return ret;
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


bool FractureUserObject::isInsideRegion2D(RealVectorValue const & point, int const i, Real & bound) const
{
    RealVectorValue localbounds;
    for (int k=0; k<_dim; ++k)
        localbounds(k)=std::max(_dimension[i](k)/2.0,bound);

    Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
    if (temp1<localbounds(0))
    {
        Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
        if (temp2<localbounds(1))
        {
            return true;
        }
    }
    return false;
}

bool FractureUserObject::isInsideRegion3D(RealVectorValue const & point, int const i, Real & bound) const
{
    RealVectorValue localbounds;
    for (int k=0; k<_dim; ++k)
        localbounds(k)=std::max(_dimension[i](k)/2.0,bound);

    //std::cout<<"I am here"<<std::endl;

    Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
    if (temp1<localbounds(0))
    {
        Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
        if (temp2<localbounds(1))
        {
            Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
            if (temp3<localbounds(2))
            {
                return true;
            }
            
        }
    }
    return false;
}


// THESE HAVE BEEN ADDED TO REFINE THE BOUNDARIES
bool FractureUserObject::isOnBoundary(Elem & elem) const
{
    for (int i=0; i<_fn; ++i)
    {
        if (isOnBoundaryOfRegion(elem,i))
            return true;
    }
    return false;
}

bool FractureUserObject::isOnBoundaryOfRegion(Elem & elem, int i) const
{
    bool ret;
    if (_dim==2)
    {
        ret=isOnBoundaryOfRegion2D(elem,i);
    }
    if (_dim==3)
    {
        ret=isOnBoundaryOfRegion3D(elem,i);
    }
    return ret;
}

bool FractureUserObject::isOnBoundaryOfRegion2D(Elem & elem, int i) const
{
    if ( elem.n_vertices()!=elem.n_nodes () )
    {
        _console<<"elem.n_vertices()!=elem.n_nodes()"<<std::endl;
        exit(1);
    }
    RealVectorValue cmax(-1e9,-1e9,-1e9);
    RealVectorValue cmin( 1e9, 1e9, 1e9);

    for (unsigned int v=0;v<elem.n_vertices(); ++v)
    {
        Point & point=elem.point(v);
        RealVectorValue rvv(point);
        for (int j=0; j<3; ++j)
        {
            cmax(j)=std::max( cmax(j),rvv(j) );
            cmin(j)=std::min( cmin(j),rvv(j) );
        }
    }

    RealVectorValue const & p0=_vertex.at(i).at(0);
    RealVectorValue const & p1=_vertex.at(i).at(1);
    RealVectorValue const & p2=_vertex.at(i).at(2);
    RealVectorValue const & p3=_vertex.at(i).at(3);

   if ( doesEdgeIntersectElement_2D(p0, p1, cmin, cmax) )
       return true;
   if ( doesEdgeIntersectElement_2D(p0, p2, cmin, cmax) )
       return true;
    if ( doesEdgeIntersectElement_2D(p1, p3, cmin, cmax) )
        return true;
     if ( doesEdgeIntersectElement_2D(p2, p3, cmin, cmax) )
         return true;

    return false;
}

bool FractureUserObject::isOnBoundaryOfRegion3D(Elem & elem, int i) const
{   
    // RealVectorValue localbounds;
    // for (int k=0; k<_dim; ++k)
    //     localbounds(k)=std::max(_dimension[i](k)/2.0,bound);

    // {
    // Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
    // if (temp1<localbounds(0))
    // {
    //     Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
    //     if (temp2<localbounds(1))
    //     {
    //         Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
    //         if ( std::fabs(temp3-localbounds(2)) <  bound+myeps )
    //         {
    //             return true;
    //         }
            
    //     }
    // }
    // }
    // {
    // Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
    // if (temp1<localbounds(0))
    // {
    //     Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
    //     if ( std::fabs( temp2- localbounds(1)) < bound+myeps )
    //     {
    //         Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
    //         if (temp3<localbounds(2))
    //         {
    //             return true;
    //         }
            
    //     }
    // }
    // }
    // {
    // Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
    // if ( std::fabs(temp1-localbounds(0)) < bound+myeps )
    // {
    //     Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
    //     if (temp2<localbounds(1))
    //     {
    //         Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
    //         if (temp3<localbounds(2))
    //         {
    //             return true;
    //         }
            
    //     }
    // }
    // }
    return false;
}
