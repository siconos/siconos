/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <iostream>

#include <string.h>
using namespace std;


CppUnit::Test* GetTest(CppUnit::TestSuite *testSuite, const std::string& seekedname)
{
  std::vector<CppUnit::Test*>::iterator testi;
  std::vector<CppUnit::Test*> allTests = testSuite->getTests();
  CppUnit::Test* returnTest;

  for (testi = allTests.begin();
       testi != allTests.end(); ++testi)
  {
    if ((*testi)->getName() == seekedname)
    {
      returnTest = *testi;
      break;
    }
  }
  return(returnTest);

};

/*! Dump a test */
int CdashDumpTest(CppUnit::Test *test, char* myname)
{
  std::cout << "ADD_TEST(" << test->getName() << " " << myname << " " << test->getName() << ")" << std::endl;
};

int CdashDump(CppUnit::Test *tests, char* myname)
{

  int count = 0;
  // try to see if this is a TestSuite
  CppUnit::TestSuite* testSuite = dynamic_cast<CppUnit::TestSuite *>(tests);
  // it's a suite, check all components
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

int ListAllTestsInTestSuite(CppUnit::Test* tests)
{
  int count = 0;
  // try to see if this is a TestSuite
  CppUnit::TestSuite* testSuite = dynamic_cast<CppUnit::TestSuite *>(tests);
  // it's a suite, check all components
  if (testSuite != NULL)
  {
    std::vector<CppUnit::Test*> allTestsVector = testSuite->getTests();
    std::vector<CppUnit::Test*>::iterator testIterator;
    for (testIterator = allTestsVector.begin();
         testIterator != allTestsVector.end();
         testIterator++)
    {
      count += ListAllTestsInTestSuite(*testIterator);
    }
  }
  else
  {
    // it's a test, get the name
    count++;
    std::cout << tests->getName() << std::endl;
  }
  return count;
}

int main(int argc, char** argv)
{

  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  CppUnit::TestSuite *testSuite = static_cast<CppUnit::TestSuite*>(registry.makeTest());

  if (argc == 2)
  {
    std::string arg = argv[1];
    if (strcmp(argv[1], "--cdash-prepare"))
    {
      CdashDump(testSuite, argv[0]);
    }
    else
    {

      // The object to run tests
      CppUnit::TextUi::TestRunner runner;

      // get the test
      CppUnit::Test * test = GetTest(testSuite, argv[1]);

      runner.addTest(test);

      bool wasSucessful = false;

      wasSucessful = runner.run("", false, true, false);

      return wasSucessful ? 0 : 1;

    }

  }
}
