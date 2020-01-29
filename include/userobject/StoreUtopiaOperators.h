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


// MOOSE includes
#include "GeneralUserObject.h"

// utopia includes
#include "utopia.hpp"
#include "utopia_fe.hpp"
#include <memory>
// Forward Declarations
class StoreUtopiaOperators;

template<>
InputParameters validParams<StoreUtopiaOperators>();

class StoreUtopiaOperators :
  public GeneralUserObject
{
public:
  typedef utopia::USparseMatrix SparseMatT;  

  // constructor
  StoreUtopiaOperators(const InputParameters & parameters);

  // returns pointer to stored transfer operators
  std::shared_ptr<SparseMatT> &
  getOperator()
  {
    return _B;
  };
 

  std::shared_ptr<SparseMatT> &
  setOperator()
  {
    return _B;
  };
    
  std::shared_ptr<void> & getVoidPointer()
    {
        return void_ptr;
    }
    
    std::shared_ptr<void> void_ptr;

  void execute(){};
  void initialize(){};
  void finalize(){};

protected:
  std::shared_ptr<SparseMatT> _B;
};


