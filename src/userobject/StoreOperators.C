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



#include "StoreOperators.h"

#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/sparse_matrix.h"

using namespace std;

registerMooseObject("parrotApp", StoreOperators);

template <>
InputParameters
validParams<StoreOperators>()
{
  InputParameters params = validParams<GeneralUserObject>();
  return params;
}

StoreOperators::StoreOperators(const InputParameters & parameters) :
    GeneralUserObject(parameters)
{
	auto &comm = _fe_problem.es().get_mesh().comm();

    // _mass_matrix = new PetscMatrix<Number>(comm);
    // _poro_mass_matrix = new PetscMatrix<Number>(comm);
    // _lump_mass_matrix = new PetscMatrix<Number>(comm);
    _interpolator = std::make_shared<PetscMatrix<Number>>(comm);
	_mass_matrix = std::make_shared<PetscMatrix<Number>>(comm);
    _poro_mass_matrix = std::make_shared<PetscMatrix<Number>>(comm);
    _lump_mass_matrix = std::make_shared<PetscMatrix<Number>>(comm);
    _poro_lump_mass_matrix = std::make_shared<PetscMatrix<Number>>(comm);
    _jac_matrix = std::make_shared<PetscMatrix<Number>>(comm);
    _bc_vec = std::make_shared<PetscVector<Number>>(comm);
    _solution_vec = std::make_shared<PetscVector<Number>>(comm);
    _value_bc_vec = std::make_shared<PetscVector<Number>>(comm);
}
