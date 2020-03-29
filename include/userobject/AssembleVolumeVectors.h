//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

#include "libmesh/petsc_vector.h"

// Forward declarations
class AssembleVolumeVectors;

template <>
InputParameters validParams<AssembleVolumeVectors>();

class AssembleVolumeVectors : public GeneralUserObject
{
public:
 AssembleVolumeVectors(const InputParameters & params);

  virtual void initialize() override {};
  virtual void execute() override ;
  virtual void finalize() override {};
  virtual void threadJoin(const UserObject & ) override {};

  std::vector< PetscVector<Number> * > const & getVolumeVectors() const {return volumes;};

  ~AssembleVolumeVectors();

  protected:

    std::vector< PetscVector<Number> * > volumes;

    bool const _hasMeshModifier;
    std::string _meshModifierName;
    bool _r_fracture;
    int _n_regions;

    unsigned _fn;

    bool _vectorAllocated;

    std::vector<int> _block_id;
    // std::vector<Real> _vector_value;

    // UserObjectName const userObjectName;
    // StoreOperators * _storeOperatorsUO;

    // std::shared_ptr<PetscMatrix<Number>> _interpolator;
    // std::shared_ptr<PetscMatrix<Number>> _mass_matrix;
    // std::shared_ptr<PetscMatrix<Number>> _lump_mass_matrix;
    // std::shared_ptr<PetscMatrix<Number>> _poro_mass_matrix;
    // std::shared_ptr<PetscMatrix<Number>> _poro_lump_mass_matrix;

    

    // void assemble_mass_matrix();
    //Real ComputeMaterialProprties(const Elem *elem);

    // bool const _constrainMatrices;
    // bool const _code_dof_map;

    // bool _hasMeshModifier;
    

    // QBase const * const & _qrule;
};

