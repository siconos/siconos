/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "ModelTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ModelTest);


void ModelTest::setUp()
{
  t0 = 0.1;
  T = 10.0;
}

void ModelTest::tearDown()
{}

void ModelTest::testBuildModel0()
{
  std::cout << "=============================" <<std::endl;
  std::cout << "=== Model tests start ...=== " <<std::endl;
  std::cout << "=============================" <<std::endl;
  std::cout << "--> Test: constructor 0." <<std::endl;
  SP::Model M(new Model(t0));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->t0() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->finalT() == -1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->currentTime() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", !M->simulation(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", !M->nonSmoothDynamicalSystem(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->title() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->author() == "nobody", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->description() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->date() == "none", true);
  std::cout << "--> Constructor 0 test ended with success." <<std::endl;
}

void ModelTest::testBuildModel1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::Model M(new Model(t0, T));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->t0() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->finalT() == T, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->currentTime() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->simulation(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->nonSmoothDynamicalSystem(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->title() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->author() == "nobody", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->description() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->date() == "none", true);
  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}

void ModelTest::testBuildModel2()
{
  std::cout << "--> Test: constructor 2." <<std::endl;
  SP::Model M(new Model(t0, T, "myModel", "SiconosTeam", "Description", "Today"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->t0() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->finalT() == T, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->currentTime() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->simulation(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->nonSmoothDynamicalSystem(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->title() == "myModel", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->author() == "SiconosTeam", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->description() == "Description", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->date() == "Today", true);
  std::cout << "--> Constructor 2 test ended with success." <<std::endl;
}
