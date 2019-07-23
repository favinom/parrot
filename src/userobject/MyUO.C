//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MyUO.h"

#include "FEProblem.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

registerMooseObject("parrotApp", MyUO);

template <>
InputParameters
validParams<MyUO>()
{
  InputParameters params = validParams<GeneralUserObject>();

  return params;
}

MyUO::MyUO(const InputParameters & parameters) :
GeneralUserObject(parameters),
_fe_problem(parameters.get<FEProblem *>("_fe_problem")),
_equationSystems(_fe_problem[0].es())
{
    std::cout<<"costruttore\n";
    
}

 void MyUO::execute()
{
    std::cout<<"execute\n";
    
    int nES=_equationSystems.n_systems();
    
    std::cout<<"The EQ has "<<nES<<std::endl;
    
    LinearImplicitSystem & lis0=_equationSystems.get_system<LinearImplicitSystem> (0);

    LinearImplicitSystem & lis1=_equationSystems.get_system<LinearImplicitSystem> (1);

    std::cout<<"System lis0 has "<<lis0.n_matrices ()<<" matrices\n";
    
    auto &_m_sys = _fe_problem[0].es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    //TransientNonlinearImplicitSystem &pippo = dynamic_cast<TransientNonlinearImplicitSystem &>(lis0);
    std::cout<<_fe_problem->getCurrentExecuteOnFlag()<<std::endl;

    PetscMatrix<Number> * petsc_mat_m =
    dynamic_cast<PetscMatrix<Number> * >(_m_sys.matrix);
//    if(_fe_problem->getCurrentExecuteOnFlag()==EXEC_LINEAR){
//        petsc_mat_m[0].set(0,0,1e15);
//        petsc_mat_m[0].close();
//    }
    std::cout<<"before writing\n";
    petsc_mat_m[0].print_matlab("daUO.m");
    std::cout<<"after writing\n";
    

    
    exit(1);
    
};

 void MyUO::initialize()
{
    std::cout<<"initialize\n";
    
}

 void MyUO::finalize()
{
    std::cout<<"finalize\n";
}

 void MyUO::threadJoin(const UserObject & uo)
{};
