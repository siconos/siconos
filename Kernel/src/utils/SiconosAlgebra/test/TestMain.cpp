#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
using namespace std;

int main()
{
  // The object to run tests
  CppUnit::TextUi::TestRunner runner;

  // Get test classes that have been registered
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();

  // Put tests into the runner
  runner.addTest(registry.makeTest());

  // Run tests
  bool wasSucessful = runner.run("", false);
  return wasSucessful;
}

