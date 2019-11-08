///* This file is part of the MOOSE framework
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

// Forward declarations
class MyUO;

template <>
InputParameters validParams<MyUO>();

/**
 * Base class for user-specific data

 */

class StoreOperators : public GeneralUserObject
{

    
    
public:
  StoreOperators(const InputParameters & params);

  virtual void execute() override ;

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize() override ;

  /**
   * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where
   * you want to do MPI communication!
   */
  virtual void finalize() override ;

  /**
   * Must override.
   *
   * @param uo The UserObject to be joined into _this_ object.  Take the data from the uo object and
   * "add" it into the data for this object.
   */
  virtual void threadJoin(const UserObject & uo) override;
 

  std::shared_ptr<PetscMatrix<Number>> &
  MassMatrix()
  {
    //std::cout << "[MASSMATRIX-StoreOperators::MassMatrix()]" << _mass_matrix.get() <<  ") (" << this << std::endl;
    return _mass_matrix;
  };


  std::shared_ptr<PetscMatrix<Number>> &
  PoroMassMatrix()
  {
    return _poro_mass_matrix;
  };

  std::shared_ptr<PetscMatrix<Number>> &
  LumpMassMatrix()
  {
    return _lump_mass_matrix;
  };


  std::shared_ptr<PetscMatrix<Number>> &
  StabMatrix()
  {
    return _stab_matrix;
  };
  
  std::shared_ptr<PetscMatrix<Number>> &
  JacMatrix()
  {
    return _jac_matrix;
  };

  std::shared_ptr<PetscMatrix<Number>> &
  PoroLumpMassMatrix()
  {
    return _poro_lump_mass_matrix;
  };



  // std::shared_ptr<PetscMatrix<Number>> &
  // setStabMatrix()
  // {
  //   return _stab_matrix;
  // };


  // std::shared_ptr<PetscMatrix<Number>> &
  // setPoroMassMatrix()
  // {
  //   return _poro_mass_matrix;
  // };
    
  // std::shared_ptr<PetscMatrix<Number>> &
  // setLumpMassMatrix()
  // {
  //   return _lump_mass_matrix;
  // };

  // std::shared_ptr<PetscMatrix<Number>> &
  // setJacMatrix()
  // {
  //   return _jac_matrix;
  // };


//   void execute(){};
//   void initialize(){};
//   void finalize(){};

protected:
  std::shared_ptr<PetscMatrix<Number>> _mass_matrix;
  std::shared_ptr<PetscMatrix<Number>> _stab_matrix;
  std::shared_ptr<PetscMatrix<Number>> _poro_mass_matrix;
  std::shared_ptr<PetscMatrix<Number>> _lump_mass_matrix;
  std::shared_ptr<PetscMatrix<Number>> _poro_lump_mass_matrix;
  std::shared_ptr<PetscMatrix<Number>> _jac_matrix;

   
    

};




