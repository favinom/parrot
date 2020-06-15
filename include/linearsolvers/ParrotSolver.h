#pragma once

#include "petscksp.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

using namespace libMesh;

class ParrotSolver
{
private:
	KSP _ksp;
	PC  _pc;

	NumericVector<Number> * _sol_NV;
 	NumericVector<Number> * _rhs_NV;
 	SparseMatrix <Number> * _mat_SM;

	PetscVector<Number> * _sol_PV;
 	PetscVector<Number> * _rhs_PV;
 	PetscMatrix<Number> * _mat_PM;

 	Vec _residual;
 	Mat _preconditioner;

 	// 1 LU
 	// 2 CHOL
 	// 3 AMG
 	int const _solverType;

 	MPI_Comm const _communicator;
 	PetscErrorCode _ierr;

 	bool _constantMatrix;
 	bool _factorized;

 	PetscErrorCode Construct();

 public:

 	ParrotSolver(int solverType);
 	ParrotSolver(int solverType, Parallel::Communicator const & _pp_comm);
 	ParrotSolver(int solverType, MPI_Comm const & comm);

 	PetscErrorCode setMatrixAndVectors(SparseMatrix <Number> * & mat_SM, NumericVector<Number> * & rhs_NV, NumericVector<Number> * & sol_NV);

 	PetscErrorCode solve();

 	PetscErrorCode computeResidual();

 	void setConstantMatrix(bool cm) {_constantMatrix=cm;}

 	~ParrotSolver();


};
