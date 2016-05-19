/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <iostream>

#include <cstring>

CppUnit::Test* GetTest(CppUnit::Test* tests, const std::string& name);
void CdashDumpTest(CppUnit::Test *test, char* myname);
int CdashDump(CppUnit::Test *tests, char* myname);

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
}

/* Dump a unit test as a cmake test */
void CdashDumpTest(CppUnit::Test *test, char* myname)
{

   std::cout << "MESSAGE( STATUS Adding unit test : " << test->getName() << " ) "
            << std::endl;
   std::cout << "ADD_CPPUNIT_TEST(" << test->getName() << " "
            << EMULATOR << " " << myname << WRAPPER << " " << test->getName()
            << ")" << std::endl;
}

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

  return count; // A verifier

}


/* <test executable> --cdash-prepare */
/* <test executable> <unit test name> */

int main(int argc, char** argv)
{

  // Registry and test suite
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  CppUnit::TestSuite *testSuite = static_cast<CppUnit::TestSuite*>(registry.makeTest());

  if (argc == 3)
  {
    std::string arg = argv[1];
    if (strcmp(argv[1], "--cdash-prepare") == 0)
    {
       std::cout << "# this is a ctest input file" << std::endl;
       std::cout << "include(SiconosTestConfig.cmake)" << std::endl;

      CdashDump(testSuite, argv[2]);
    }
  }
  else if (argc == 2)
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

  else
  {
    std::cerr << "Error, no test given" << std::endl;
    return 1;
  }
}
