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

#include "Interaction.hpp"
#include "FirstOrderType1RTest.hpp"
#include "SSLH.hpp"



#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);


// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderType1RTest);


void FirstOrderType1RTest::setUp()
{
}

void FirstOrderType1RTest::tearDown()
{
}

// data constructor
void FirstOrderType1RTest::testBuildFirstOrderType1R1()
{
  std::cout << "======================================" <<std::endl;
  std::cout << "=== FirstOrderType1R tests start ...== " <<std::endl;
  std::cout << "======================================" <<std::endl;
  std::cout << "--> Test: constructor 1 " <<std::endl;

  SP::FirstOrderType1R R1(new FirstOrderType1R("TestPlugin:hT1", "TestPlugin:gT1"));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1b : ", R1->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1c : ", R1->getSubType() == RELATION::Type1R, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1d : ", R1->gethName()=="TestPlugin:hT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1e : ", R1->getgName()=="TestPlugin:gT1", true);
  std::string plugin = "TestPlugin" + SSLH::getSharedLibraryExtension();
  R1->setComputeJachxFunction(plugin, "Jh0T1");
  R1->setComputeJacglambdaFunction(plugin, "Jg0T1");
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1e : ", R1->getJachName(0)=="TestPlugin:Jh0T1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1g : ", R1->getJacgName(0)=="TestPlugin:Jg0T1", true);
  std::cout << "--> Constructor1 test ended with success." <<std::endl;
}

void FirstOrderType1RTest::testBuildFirstOrderType1R2()
{
  std::cout << "--> Test: constructor data (2)." <<std::endl;
  SP::FirstOrderType1R R2(new FirstOrderType1R("TestPlugin:hT1", "TestPlugin:gT1", "TestPlugin:Jh0T1", "TestPlugin:Jg0T1"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2b : ", R2->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2c : ", R2->getSubType() == RELATION::Type1R, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2d : ", R2->gethName()=="TestPlugin:hT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2e : ", R2->getgName()=="TestPlugin:gT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2e : ", R2->getJachName(0)=="TestPlugin:Jh0T1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2g : ", R2->getJacgName(0)=="TestPlugin:Jg0T1", true);
  std::cout << "--> Constructor2 test ended with success." <<std::endl;
}
