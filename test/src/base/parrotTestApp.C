//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "parrotTestApp.h"
#include "parrotApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<parrotTestApp>()
{
  InputParameters params = validParams<parrotApp>();
  return params;
}

parrotTestApp::parrotTestApp(InputParameters parameters) : MooseApp(parameters)
{
  parrotTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

parrotTestApp::~parrotTestApp() {}

void
parrotTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  parrotApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"parrotTestApp"});
    Registry::registerActionsTo(af, {"parrotTestApp"});
  }
}

void
parrotTestApp::registerApps()
{
  registerApp(parrotApp);
  registerApp(parrotTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
parrotTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  parrotTestApp::registerAll(f, af, s);
}
extern "C" void
parrotTestApp__registerApps()
{
  parrotTestApp::registerApps();
}
