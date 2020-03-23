#ifndef _KSPPARROTIMPL_H
#define _KSPPARROTIMPL_H

#include "FEProblem.h"
#include <petsc/private/kspimpl.h>
#include "libmesh/petsc_vector.h"
//#include <petsc/private/snesimpl.h>

typedef struct {
    PC * local_pc;
    int * factorized;
    FEProblem * fe_problem;
    //PetscVector<Number> &sol_vec;
} KSP_PARROT;

#endif
