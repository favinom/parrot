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

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#include "StoreUtopiaOperators.h"
#include "MooseVariable.h"

registerMooseObject("parrotApp", StoreUtopiaOperators);

template<>
InputParameters validParams<StoreUtopiaOperators>()
{
  InputParameters params = validParams<GeneralUserObject>();
  return params;
}

StoreUtopiaOperators::StoreUtopiaOperators(const InputParameters & parameters) :
    GeneralUserObject(parameters)
{
   _B           = std::make_shared<SparseMatT>();
}
