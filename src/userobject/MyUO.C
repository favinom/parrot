//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MyUO.h"

#include "libmesh/sparse_matrix.h"

registerMooseObject("parrotApp", MyUO);

template <>
InputParameters
validParams<MyUO>()
{
  InputParameters params = validParams<GeneralUserObject>();

  return params;
}

MyUO::MyUO(const InputParameters & parameters) :
GeneralUserObject(parameters)
{
    std::cout<<"costruttore\n";
    
}

 void MyUO::execute()
{
    std::cout<<"execute\n";
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
