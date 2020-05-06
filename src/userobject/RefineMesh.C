//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RefineMesh.h"

// #include "FEProblem.h"
// #include "libmesh/nonlinear_implicit_system.h"
// #include "libmesh/petsc_matrix.h"
// #include "libmesh/sparse_matrix.h"
// #include "libmesh/equation_systems.h"
// #include "libmesh/linear_implicit_system.h"
// #include "libmesh/transient_system.h"

registerMooseObject("parrotApp", RefineMesh);

template <>
InputParameters
validParams<RefineMesh>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  params.addRequiredParam< std::vector<int> >("refinements","refinements");
  params.addRequiredParam< std::string >("filename","filename");
  return params;
}

RefineMesh::RefineMesh(const InputParameters & parameters) :
GeneralUserObject(parameters),
_fe_problem(parameters.get<FEProblem *>("_fe_problem")),
_equationSystems(_fe_problem->es()),
_mesh(_fe_problem->mesh()),
_meshBase( _equationSystems.get_mesh() ),
_pp_comm(_mesh.comm()),
_hasMeshModifier( isParamValid("fractureMeshModifier") ),
_refinements( parameters.get< std::vector<int> >("refinements") ),
_filename( parameters.get< std::string >("filename") )
//
{
    if (_hasMeshModifier)
        _meshModifierName=getParam<std::string>("fractureMeshModifier");

    if ( _refinements.size()%2 != 0 )
    {
        std::cout<<"The size of refinements has to be even, exiting...\n";
        exit(1);
    }
    else
    {
        std::cout<<"The size of refinements is "<<_refinements.size()<<"\n";
    }
    std::cout<<"We are going to perform this sequence of refinements\n";
    for (int i=0; i<_refinements.size(); ++++i)
    {
        _console<<_refinements.at(i)<<" AR\n";
        _console<<_refinements.at(i+1)<<" UR\n";
    }
}

void RefineMesh::initialize()
{
    //DistributedMesh mesh(_pp_comm);
    //MeshTools::Generation::build_cube(mesh,9,20,9,0,1,0,2.25,0,1,HEX8);
    //MeshTools::Generation::build_square (mesh, 5, 5);

    DistributedMesh & distributedMesh( dynamic_cast<DistributedMesh &>(_meshBase) );
    _distributedMesh = &distributedMesh;
    _meshRefinement = new MeshRefinement(_distributedMesh[0]);

    for (int i=0 ;i<_refinements.size(); ++++i)
    {
        std::cout<<"\nAdaptivty loop " <<i/2<<": "<<_refinements.at(i)<<" AMR and "<<_refinements.at(i+1)<<" UMR"<<std::endl;

        for (int j=0; j<_refinements.at(i); ++j)
        {
            std::cout<<" Doing AMR "<<j+1<<std::endl;
            doAMR();
        }
        std::cout<<" Doing "<<_refinements.at(i+1)<<" UMR"<<std::endl;
        if (_refinements.at(i+1)>0)
            doUMR(_refinements.at(i+1));
    }

    std::cout<<"\nWriting mesh file\n";
    _distributedMesh[0].write(_filename);
    std::cout<<"Done!\n";
}

  void RefineMesh::doAMR()
  {
    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );

    MeshBase::const_element_iterator     el=_distributedMesh->active_local_elements_begin();
    MeshBase::const_element_iterator end_el=_distributedMesh->active_local_elements_end  ();

    std::cout<<"  Looping over elements\n";
    for ( ; el != end_el ; ++el)
    {
        Elem * elem = *el;
        Point center=elem->centroid();
        Real hmax=elem->hmax();
        hmax=hmax/2.0*std::sqrt(2.0);
        if (_fractureUserObject.isInside(center,hmax))
            elem[0].set_refinement_flag(Elem::REFINE);
    }

    std::cout<<"  Actual refinement\n";
    _meshRefinement->refine_elements();
    std::cout<<"  Done!\n";
  }

  void RefineMesh::doUMR(int i)
  {
    _meshRefinement->uniformly_refine(i);
  }

//refine_and_coarsen_elements();

// mesh.contract();
// mesh.prepare_for_use();

// mesh.redistribute ();
// mesh.update_post_partitioning ();
