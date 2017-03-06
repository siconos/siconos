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
#include "NonSmoothDynamicalSystemTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NonSmoothDynamicalSystemTest);


void NonSmoothDynamicalSystemTest::setUp()
{


void NonSmoothDynamicalSystemTest::tearDown()
{}

// copy constructor
void NonSmoothDynamicalSystemTest::testBuildNonSmoothDynamicalSystem2()
{
  std::cout << "------- Copy Constructor test -------" <<std::endl;
  SP::NonSmoothDynamicalSystem  nsds1(new NonSmoothDynamicalSystem(tmpxml));
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(*nsds1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2A : ", nsds->getDSVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2B : ", nsds->dynamicalSystem(0)->number() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2C : ", nsds->dynamicalSystem(1)->number() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2D : ", nsds->getInteractionVectorSize() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->interaction(0)->number() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2F : ", nsds->isBVP() == false, true);
  std::cout << "------- Constructor copy NonSmoothDynamicalSystem ok -------" <<std::endl;
}


// insertDynamicalSystem
void NonSmoothDynamicalSystemTest::testinsertDynamicalSystem()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(tmpxml));
  xmlNode *node2 = SiconosDOMTreeTools::findNodeChild(node, "DS_Definition");
  xmlNode * node3 = SiconosDOMTreeTools::findNodeChild(node2, "LagrangianLinearTIDS");
  SP::DynamicalSystemXML tmpdsxml(new LagrangianLinearTIDSXML(node3, false));

  SP::DynamicalSystem ltids(new LagrangianLinearTIDS(tmpdsxml));
  ltids ->setNumber(23);

  nsds->insertDynamicalSystem(ltids);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testinsertDynamicalSystemA : ", nsds->getDSVectorSize() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemB: ", nsds->dynamicalSystem(0)->number() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testinsertDynamicalSystemC : ", nsds->dynamicalSystem(1)->number() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testinsertDynamicalSystemC : ", nsds->dynamicalSystem(2)->number() == 23, true);
  std::cout << "------- test insertDynamicalSystem ok -------" <<std::endl;
}

// insertInteraction
void NonSmoothDynamicalSystemTest::testinsertInteraction()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(tmpxml));
  xmlNode *node2 = SiconosDOMTreeTools::findNodeChild(node, "Interaction_Definition");
  xmlNode *node3 = SiconosDOMTreeTools::findNodeChild(node2, "Interaction");
  vector<int> tmp;
  tmp.resize(2, 1);
  SP::InteractionXML interxml(new InteractionXML(node3, tmp));
  SP::Interaction inter(new Interaction(interxml));
  nsds->insertInteraction(inter);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2D : ", nsds->getInteractionVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->interaction(0)->number() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->interaction(1)->number() == 12, true);
  std::cout << " ------- test insertInteractiontest ok -------" <<std::endl;
}

