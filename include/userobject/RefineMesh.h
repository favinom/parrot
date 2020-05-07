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
  MooseMesh & _mooseMesh;
  Parallel::Communicator const & _pp_comm;
  MeshBase  & _meshBase;
  UnstructuredMesh * _unstructuredMesh;
  DistributedMesh  * _distributedMesh;
  MeshRefinement   * _meshRefinement;

  Mesh * _mesh;

  bool _hasMeshModifier;
  int _flag;
  std::string _meshModifierName;
  std::vector<int> const _refinements;
  std::string  const     _filename;

    
    
public:
  RefineMesh(const InputParameters & params);

  virtual void execute() override {};
  virtual void finalize() override {};
  virtual void threadJoin(const UserObject &) override {};

  virtual void initialize() override;
  void doAMR();
  void getInfo(MeshBase * mesh);
  void doUMR(int i);

  void doRefine(UnstructuredMesh & mesh, std::vector<int> const & refinements);

};

