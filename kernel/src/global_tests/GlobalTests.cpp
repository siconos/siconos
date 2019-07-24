/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "GlobalTests.hpp"



GlobalTest::GlobalTest(const std::string& logName): numberOfFailedTests(0)
{
  out.open(logName.c_str());
}

GlobalTest::~GlobalTest()
{
  out.close();
  std::set<FunctionTest*>::iterator it;
  for (it = testsList.begin(); it != testsList.end(); ++it)
    if ((*it)) delete(*it);
  testsList.clear();
}

void GlobalTest::addTest(const std::string& name, func f)
{
  std::cout << "Add a new test: " << name <<std::endl <<std::endl;
  testsList.insert(new FunctionTest(name, f));
}
void GlobalTest::run()
{
  std::set<FunctionTest*>::iterator it;
  for (it = testsList.begin(); it != testsList.end(); ++it)
  {
    (*it)->run();
    if ((*it)->hasFailed()) numberOfFailedTests++;
  }
}

void GlobalTest::print()
{
  out << " ======== This is a log file that contains global tests results and comments ======== " <<std::endl;

  out << "Number of test function run: " << testsList.size() <<std::endl;
  out << "Number of failed tests: " << numberOfFailedTests <<std::endl;

  out << "The following tests failed: " <<std::endl;
  std::set<FunctionTest*>::iterator it;
  for (it = testsList.begin(); it != testsList.end(); ++it)
  {
    if ((*it)->hasFailed()) out  << "  - " << (*it)->getName() <<std::endl;
  }
  out << "======================================================================================= " <<std::endl;

}

