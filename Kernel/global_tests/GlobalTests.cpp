/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
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

#include "GlobalTests.h"

using namespace std;

GlobalTest::GlobalTest(const std::string& logName): numberOfFailedTests(0)
{
  out.open(logName.c_str());
}

GlobalTest::~GlobalTest()
{
  out.close();
  set<FunctionTest*>::iterator it;
  for (it = testsList.begin(); it != testsList.end(); ++it)
    if ((*it) != NULL) delete(*it);
  testsList.clear();
}

void GlobalTest::addTest(const string& name, func f)
{
  cout << "Add a new test: " << name << endl << endl;
  testsList.insert(new FunctionTest(name, f));
}
void GlobalTest::run()
{
  set<FunctionTest*>::iterator it;
  for (it = testsList.begin(); it != testsList.end(); ++it)
  {
    (*it)->run();
    if ((*it)->hasFailed()) numberOfFailedTests++;
  }
}

void GlobalTest::print()
{
  out << " ======== This is a log file that contains global tests results and comments ======== " << endl;

  out << "Number of test function run: " << testsList.size() << endl;
  out << "Number of failed tests: " << numberOfFailedTests << endl;

  out << "The following tests failed: " << endl;
  set<FunctionTest*>::iterator it;
  for (it = testsList.begin(); it != testsList.end(); ++it)
  {
    if ((*it)->hasFailed()) out  << "  - " << (*it)->getName() << endl;
  }
  out << "======================================================================================= " << endl;

}

