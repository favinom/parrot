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

#include "libmesh/mesh_tools.h"

registerMooseObject("parrotApp", RefineMesh);

template <>
InputParameters
validParams<RefineMesh>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  params.addRequiredParam< std::vector<int> >("refinements","refinements");
  params.addRequiredParam< std::string >("filename","filename");
  params.addRequiredParam< int >("flag","flag");
  return params;
}

RefineMesh::RefineMesh(const InputParameters & parameters) :
GeneralUserObject(parameters),
_fe_problem(parameters.get<FEProblem *>("_fe_problem")),
_equationSystems(_fe_problem->es()),
_mooseMesh(_fe_problem->mesh()),
_pp_comm(_mooseMesh.comm()),
_meshBase( _equationSystems.get_mesh() ),
_hasMeshModifier( isParamValid("fractureMeshModifier") ),
_flag( parameters.get< int >("flag") ),
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
}

void RefineMesh::initialize()
{
    _unstructuredMesh=dynamic_cast<UnstructuredMesh *>(&_meshBase);
    if (!_hasMeshModifier)
    {
        std::cout<<"\nWriting mesh file\n";
        _unstructuredMesh->write(_filename);
        std::cout<<"Done!\n";
        return;
    }

    if (_flag==1)
    {
        _distributedMesh=dynamic_cast<DistributedMesh *>(_unstructuredMesh);
    }
    if (_flag==2)
    {
        _distributedMesh=new DistributedMesh(_pp_comm);
        _distributedMesh->read("bcMesh.xdr");
    }
    if (_flag==3)
    {
        _distributedMesh=new DistributedMesh(_unstructuredMesh[0]);
    }
    _meshRefinement=new MeshRefinement(_distributedMesh[0]);


    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );

    MeshBase::const_element_iterator     el=_distributedMesh->active_elements_begin();
    MeshBase::const_element_iterator end_el=_distributedMesh->active_elements_end  ();
    std::cout<<"  Looping over elements\n";
    for ( ; el != end_el ; ++el)
    {
        Elem * elem = *el;
        Point center=elem->centroid();
        Real hmax=elem->hmax();
        hmax=hmax/2.0*std::sqrt(2.0);
        if (_fractureUserObject.isInside(center,hmax))
            elem[0].set_refinement_flag(Elem::REFINE);
        else
            elem[0].set_refinement_flag(Elem::DO_NOTHING);
    }

    std::cout<<"  Actual refinement\n";
    _meshRefinement->refine_elements();
    std::cout<<"  Done!\n";

    _distributedMesh->allgather();
    getInfo(_distributedMesh);


    el=_distributedMesh->active_elements_begin();
    end_el=_distributedMesh->active_elements_end  ();

    std::cout<<"  Looping over elements\n";
    for ( ; el != end_el ; ++el)
    {
        Elem * elem = *el;
        Point center=elem->centroid();
        Real hmax=elem->hmax();
        hmax=hmax/2.0*std::sqrt(2.0);
        if (_fractureUserObject.isInside(center,hmax))
            elem[0].set_refinement_flag(Elem::REFINE);
        else
            elem[0].set_refinement_flag(Elem::DO_NOTHING);
    }

    std::cout<<"  Actual refinement\n";
    _meshRefinement->refine_elements();
    std::cout<<"  Done!\n";

    delete _meshRefinement;
    if (_flag==2 || _flag==3)
    {
        delete _distributedMesh;
    }
}

  void RefineMesh::doAMR()
  {
    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );

    MeshBase::const_element_iterator     el=_mesh->elements_begin();
    MeshBase::const_element_iterator end_el=_mesh->elements_end  ();

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
  void RefineMesh::getInfo(MeshBase * mesh)
  {
    MeshBase::const_element_iterator     el=mesh->active_local_elements_begin();
    MeshBase::const_element_iterator end_el=mesh->active_local_elements_end  ();
    int counter = 0;
    for ( ; el != end_el ; ++el)
    {
        counter++;
    }
    std::cout<<counter<<std::endl;

        el=mesh->local_elements_begin();
    end_el=mesh->local_elements_end  ();
    
    counter = 0;
    for ( ; el != end_el ; ++el)
    {
        counter++;
    }
    std::cout<<counter<<std::endl;

    el=mesh->active_elements_begin();
    end_el=mesh->active_elements_end  ();
    counter = 0;
    for ( ; el != end_el ; ++el)
    {
        counter++;
    }
    std::cout<<counter<<std::endl;
}

void RefineMesh::doRefine(UnstructuredMesh & mesh, std::vector<int> const & refinements)
{
    std::cout<<"We are going to perform this sequence of refinements\n";
    for (int i=0; i<refinements.size(); ++++i)
    {
        std::cout<<refinements.at(i)<<" AR\n";
        std::cout<<refinements.at(i+1)<<" UR\n";
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


// //MeshTools::Generation::build_cube(mesh,9,20,9,0,1,0,2.25,0,1,HEX8);

    // //MeshTools::Generation::build_square (mesh, 5, 5);


//refine_and_coarsen_elements();

// mesh.contract();
// mesh.prepare_for_use();

// mesh.redistribute ();
// mesh.update_post_partitioning ();

    // _distributedMesh=dynamic_cast<DistributedMesh *>(&_meshBase);
    // _mesh=_distributedMesh;
    // _meshRefinement=new MeshRefinement(_mesh[0]);

    // getInfo(_mesh);

    // doRefine(_mesh[0], _refinements);

    // std::cout<<"\nWriting mesh file\n";
    // _mesh->write(_filename);
    // std::cout<<"Done!\n";

    // delete _mesh;
    // delete _meshRefinement;
    //  virtual void clear_extra_ghost_elems() { _extra_ghost_elems.clear(); }

//    libMesh::MeshTools::correct_node_proc_ids(_distributedMesh[0]);

