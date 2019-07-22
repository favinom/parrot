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
#include "MinMaxRegion.h"

registerMooseObject("parrotApp", RegionMaterial);

template<>
InputParameters validParams<RegionMaterial>()
{
    InputParameters params = validParams<Material>();
    params.addRequiredParam<std::string>("regionMeshModifier",
                                         "regionMeshModifier");
    return params;
}

RegionMaterial::RegionMaterial(const InputParameters & parameters) :
Material(parameters),
_regionID(declareProperty<int>("RegionID")),
_regionIDReal(declareProperty<Real>("RegionIDReal")),
_meshModifierName(getParam<std::string>("regionMeshModifier"))
{}

void
RegionMaterial::computeQpProperties()
{
    
    MeshModifier const & _myMeshModifier= _app.getMeshModifier(_meshModifierName.c_str());
    //MeshModifier const * _myMeshModifier_ptr= &_app.getMeshModifier(_meshModifierName.c_str());
    
    std::vector<std::string> test;
    test=_app.getMeshModifierNames();
     RegionUserObject const & regionUserObject=dynamic_cast<RegionUserObject const &>(_myMeshModifier);
    
    _regionID[_qp]=-1;
    std::vector<int> ciao;
    ciao.clear();

    ciao=regionUserObject.whichIsInside(_q_point[_qp]);
    if (ciao.size()!=1)
    {
        std::cout<<"maggiore\n";
        exit(1);
    }
    _regionID[_qp]=ciao.at(0);
    _regionIDReal[_qp]=(Real)_regionID[_qp];
}
