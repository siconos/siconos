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
#include "LagrangianRheonomousRTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianRheonomousRTest);


void LagrangianRheonomousRTest::setUp()
{}


void LagrangianRheonomousRTest::tearDown()
{}

// data constructor:
void LagrangianRheonomousRTest::testBuildLagrangianRheonomousR0()
{
  SP::LagrangianRheonomousR R1(new LagrangianRheonomousR("TestPlugin:hRheo", "TestPlugin:G0Rheo", "TestPlugin:hDot"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianRheonomousR3a : ", R1->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianRheonomousR3b : ", R1->getSubType() == RELATION::RheonomousR, true);
  std::cout << " data Constructor LagrangianRheonomousR ok" <<std::endl;
}

