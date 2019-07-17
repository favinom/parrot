#include "parrotApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<parrotApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

parrotApp::parrotApp(InputParameters parameters) : MooseApp(parameters)
{
  parrotApp::registerAll(_factory, _action_factory, _syntax);
}

parrotApp::~parrotApp() {}

void
parrotApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"parrotApp"});
  Registry::registerActionsTo(af, {"parrotApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
parrotApp::registerApps()
{
  registerApp(parrotApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
parrotApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  parrotApp::registerAll(f, af, s);
}
extern "C" void
parrotApp__registerApps()
{
  parrotApp::registerApps();
}
