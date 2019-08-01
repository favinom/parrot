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

#include "RegionMaterial.h"

#include "RegionUserObject.h"
//#include "MinMaxRegion.h"

registerMooseObject("parrotApp", RegionMaterial);

template<>
InputParameters validParams<RegionMaterial>()
{
    InputParameters params = validParams<Material>();
    params.addRequiredParam<std::string>("regionMeshModifier",
                                         "regionMeshModifier");
    
    params.addParam<bool>("countFractures","countFractures");
    
    return params;
}

RegionMaterial::RegionMaterial(const InputParameters & parameters) :
Material(parameters),
_regionID(declareProperty<int>("RegionID")),
_regionIDReal(declareProperty<Real>("RegionIDReal")),
_howManyFractures(declareProperty<Real>("HowManyFractures")),
_meshModifierName(getParam<std::string>("regionMeshModifier"))
{
    
    if (parameters.isParamValid("countFractures"))
    {
        _countFractures=getParam<bool>("countFractures");
    }
    else
        _countFractures=false;
    
    _elem_counter=0;
}

void
RegionMaterial::computeQpProperties()
{
    
    MeshModifier const & _myMeshModifier= _app.getMeshModifier(_meshModifierName.c_str());
    
//    std::vector<std::string> test;
//    test=_app.getMeshModifierNames();
     RegionUserObject const & regionUserObject=dynamic_cast<RegionUserObject const &>(_myMeshModifier);
    
    
    std::vector<int> ciao;
    ciao.clear();

    if (_countFractures)
    {
        _regionID[_qp]=_current_elem->id();
        _regionIDReal[_qp]=(Real)_regionID[_qp];
        ciao=regionUserObject.whichIsInside(_q_point[_qp]);
    }
    else
    {
        _regionID[_qp]=-1;
        ciao=regionUserObject.whichIsInside(_q_point[_qp]);
        if (ciao.size()!=1)
        {
            std::cout<<"maggiore\n";
            exit(1);
        }
        _regionID[_qp]=ciao.at(0);
        _regionIDReal[_qp]=(Real)_regionID[_qp];
    }
    _howManyFractures[_qp]=ciao.size();
    
    
    if (_qp+1 == _qrule->n_points() )
    {
        //std::cout<<"ho fatto tutto l'elemento "<<_elem_counter<<"\n";
        _elem_counter++;
        
        Real mmin=1000;
        Real mmax=-1000;
        
        for (int qp=0; qp<_qrule->n_points(); ++qp)
        {
            mmin=std::min(_howManyFractures[qp],mmin);
            mmax=std::max(_howManyFractures[qp],mmax);
        }

        
        if ( std::fabs(mmax-mmin)>0.2 )
        {
            std::cout<<"min="<<mmin<<std::endl;
            std::cout<<"max="<<mmax<<std::endl;
            std::cout<<_current_elem[0]<<std::endl;
            std::cout<<"\n\n\n";
        }
        
    }
    
}
