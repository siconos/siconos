
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

using namespace std;

#include <iostream>
int main()
{
  // on déclare un runner
  CppUnit::TextUi::TestRunner runner;

  // on récupère les classes de tests déclarées dans le registry.
  // chaque classe de test doit bien sur se referencer dans le registry
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();

  // on ajoute tous les tests dans le runner
  runner.addTest(registry.makeTest());

  // on lance les tests
  runner.run("", false, true, false);
}
