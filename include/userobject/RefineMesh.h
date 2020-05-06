//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
//#include "UserObject.h"

#include "GeneralUserObject.h"

#include "libmesh/distributed_mesh.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"

// Forward declarations
class RefineMesh;

template <>
InputParameters validParams<RefineMesh>();

/**
 * Base class for user-specific data
 */
class RefineMesh : public GeneralUserObject
{
protected:

  FEProblem * _fe_problem;
  EquationSystems & _equationSystems;
  MooseMesh & _mesh;
  MeshBase  & _meshBase;
  Parallel::Communicator const & _pp_comm;
  bool _hasMeshModifier;
  std::string _meshModifierName;
  std::vector<int> const _refinements;
  DistributedMesh * _distributedMesh;
  MeshRefinement  * _meshRefinement;
  std::string  const     _filename;

    
    
public:
  RefineMesh(const InputParameters & params);

  virtual void execute() override {};

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize() override ;

  /**
   * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where
   * you want to do MPI communication!
   */
  virtual void finalize() override {};

  /**
   * Must override.
   *
   * @param uo The UserObject to be joined into _this_ object.  Take the data from the uo object and
   * "add" it into the data for this object.
   */
  virtual void threadJoin(const UserObject &) override {};

  void doAMR();
  void doUMR(int i);

};

