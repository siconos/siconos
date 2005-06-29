#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include "SiconosException.h"
#include <iostream>
using namespace std;

int main()
{
  try
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

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught" << endl;
  }
}

