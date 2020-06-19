//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FractureRefinement.h"

registerMooseObject("MooseApp", FractureRefinement);

template <>
InputParameters
validParams<FractureRefinement>()
{
    
    InputParameters params = validParams<MeshModifier>();
    params.addClassDescription("Changes the subdomain ID of elements either (XOR) inside or outside "
                               "the specified box to the specified ID.");    
    params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");
    params.addParam< bool >("doBoundaryRefinement","doBoundaryRefinement");
    params.addParam< std::vector<int> >("refinements","refinements");
    params.addParam<std::string>("outputFileName","outputFileName");
    
    
    return params;
}

FractureRefinement::FractureRefinement(const InputParameters & parameters) :
MeshModifier(parameters),
_meshModifierName(getParam<std::string>("fractureMeshModifier")),
_myMeshModifier(_app.getMeshModifier(_meshModifierName.c_str())),
_fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) ),
_hasRefinementVector( isParamValid("refinements") ),
_hasOutputFileName( isParamValid("outputFileName") )
{
	if (_hasRefinementVector)
	{
		_refinements=getParam<std::vector<int> >("refinements");
		if ( _refinements.size()%2 != 0 )
    	{
        	std::cout<<"The size of refinements has to be even, exiting...\n";
        	exit(1);
    	}
	}
	else
		std::cout<<"Maybe print mesh\n";

	if (_hasOutputFileName)
	{
		_outputFileName=getParam<std::string >("outputFileName");
	}
  if ( isParamValid("doBoundaryRefinement") )
  {
    _doBoundary=getParam<bool>("doBoundaryRefinement");
  }
  else
    _doBoundary=false;
}

void FractureRefinement::modify()
{
	_meshBase=&_mesh_ptr->getMesh();
	//Parallel::Communicator const & _pp_comm( _mesh_ptr->getMesh().comm() );
	distributedMesh=dynamic_cast<UnstructuredMesh *>(_meshBase);
	meshRefinement=new MeshRefinement(distributedMesh[0]);

	if (_hasRefinementVector)
		doRefine(_refinements);

	if (_hasOutputFileName)
	{
		std::cout<<"\nWriting mesh file\n";
		distributedMesh->write(_outputFileName);
		std::cout<<"Done!\n";
	}
}

void FractureRefinement::doAMR()
 {
 	MeshBase::const_element_iterator     el=distributedMesh->active_elements_begin();
 	MeshBase::const_element_iterator end_el=distributedMesh->active_elements_end  ();

 	std::cout<<"  Looping over elements\n";
 	for ( ; el != end_el ; ++el)
 	{
 		Elem * elem = *el;
   		Point center=elem->centroid();
   		Real hmax=elem->hmax();
   		hmax=hmax/2.0*std::sqrt(2.0);
      if (_doBoundary)
      {
        if ( _fractureUserObject.isOnBoundary(elem[0]) )
        {
          elem[0].set_refinement_flag(Elem::REFINE);
        }
      }
      else
      {
        if (_fractureUserObject.isInside(center,0.0) || _fractureUserObject.isOnBoundary(elem[0]))
          elem[0].set_refinement_flag(Elem::REFINE);
      }
   		//else
   		//	elem[0].set_refinement_flag(Elem::DO_NOTHING);
   	}
   	std::cout<<"  Actual refinement\n";
   	meshRefinement->refine_elements();
   	std::cout<<"  Done!\n";
}

  void FractureRefinement::doUMR(int i)
  {
    meshRefinement->uniformly_refine(i);
  }


void FractureRefinement::doRefine(std::vector<int> const & refinements)
{
    std::cout<<"We are going to perform this sequence of refinements\n";
    for (int i=0; i<refinements.size(); ++++i)
    {
        std::cout<<refinements.at(i)<<" AMR and ";
        std::cout<<refinements.at(i+1)<<" UMR\n";
    }

    for (int i=0;i<refinements.size(); ++++i)
    {
        std::cout<<"\nAdaptivty loop " <<i/2<<": "<<refinements.at(i)<<" AMR and "<<refinements.at(i+1)<<" UMR"<<std::endl;
        for (int j=0; j<refinements.at(i); ++j)
        {
            std::cout<<" Doing AMR "<<j+1<<std::endl;
            doAMR();
        }

        std::cout<<" Doing "<<refinements.at(i+1)<<" UMR"<<std::endl;
        if (refinements.at(i+1)>0)
            doUMR(refinements.at(i+1));
    }

}
