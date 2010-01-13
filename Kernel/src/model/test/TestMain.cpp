/* Siconos-Kernel, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <iostream>

#include <string.h>
using namespace std;

/* get a test pointer in the test suite */

CppUnit::Test* GetTest(CppUnit::Test* tests, const std::string& name)
{

  CppUnit::TestSuite* testSuite = dynamic_cast<CppUnit::TestSuite *>(tests);

  CppUnit::Test* returnTest = NULL;

  if (testSuite)
  {
    if (testSuite->getName() == name)
    {
      return (testSuite);
    }
    else
    {
      std::vector<CppUnit::Test*> allTests = testSuite->getTests();
      std::vector<CppUnit::Test*>::iterator testi;
      for (testi = allTests.begin();
           testi != allTests.end();
           testi++)
      {
        returnTest = GetTest(*testi, name);
        if (returnTest)
          return (returnTest);
      }
    }
  }
  else
  {
    if (tests->getName() == name)
    {
      return (tests);
    }
  }
  return NULL;
};

/* Dump a unit test as a cmake test */
int CdashDumpTest(CppUnit::Test *test, char* myname)
{

  std::cout << "MESSAGE( STATUS Adding unit test : " << test->getName() << " ) "
            << std::endl;

  std::cout << "ADD_TEST(" << test->getName() << " "
            << myname << ".ldwrap " << test->getName() << ")"
            << std::endl;
};

/* Dump the test suite */
int CdashDump(CppUnit::Test *tests, char* myname)
{

  int count = 0;

  CppUnit::TestSuite* testSuite = dynamic_cast<CppUnit::TestSuite *>(tests);

  if (testSuite != NULL)
  {

    std::vector<CppUnit::Test*> allTests = testSuite->getTests();
    std::vector<CppUnit::Test*>::iterator testi;
    for (testi = allTests.begin();
         testi != allTests.end();
         testi++)
    {
      count += CdashDump(*testi, myname);
    }

  }
  else
  {
    CdashDumpTest(tests, myname);
  }



};


/* <test executable> --cdash-prepare */
/* <test executable> <unit test name> */

int main(int argc, char** argv)
{

  // Registry and test suite
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  CppUnit::TestSuite *testSuite = static_cast<CppUnit::TestSuite*>(registry.makeTest());

  if (argc == 2)
  {
    std::string arg = argv[1];
    if (strcmp(argv[1], "--cdash-prepare") == 0)
    {
      CdashDump(testSuite, argv[0]);
    }
    else
    {

      // The object to run tests
      CppUnit::TextUi::TestRunner runner;

      // get the test
      CppUnit::Test * test = GetTest(testSuite, argv[1]);

      if (test != NULL)
      {
        runner.addTest(test);

        bool wasSucessful = false;

        wasSucessful = runner.run("", false, true, false);

        return wasSucessful ? 0 : 1;
      }
      else
      {
        std::cerr << "Cannot find test : " << argv[1] << std::endl;
        return 1;
      }

    }

  }
}
