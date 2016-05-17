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
