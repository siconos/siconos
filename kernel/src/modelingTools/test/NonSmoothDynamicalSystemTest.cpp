/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "LagrangianLinearTIR.hpp"
#include "NewtonImpactNSL.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NonSmoothDynamicalSystemTest);


void NonSmoothDynamicalSystemTest::setUp()
{}


void NonSmoothDynamicalSystemTest::tearDown()
{}


// insertDynamicalSystem
void NonSmoothDynamicalSystemTest::testinsertDynamicalSystem()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(1., 10.));
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" test ndsn build: ", nsds->currentTime() == 1., true);

  SP::DynamicalSystem ds(new LagrangianDS(std::make_shared<SiconosVector>(3),
                                          std::make_shared<SiconosVector>(3)));
  ds->setNumber(23);

  try
  {
    SP::DynamicalSystem dsnull;
    nsds->insertDynamicalSystem(dsnull);
  }
  catch(const siconos::exception& e)
  {
    /*  Pass */
    std::cout << "testinsertDynamicalSystemNull: success!" << std::endl;
  }
  catch(const std::exception& e)
  {
    std::cout << "testinsertDynamicalSystemNull:" << e.what() << std::endl;
    CPPUNIT_FAIL("testinsertDynamicalSystemNull: unexpected exception ");
  }

  nsds->insertDynamicalSystem(ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemA: ", nsds->getNumberOfDS() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemC: ", nsds->dynamicalSystem(23)->number() == 23, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemD: ", nsds->topology()->hasDynamicalSystem(ds) == true, true);

  // Try again : must be ignored
  nsds->insertDynamicalSystem(ds);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemE: ", nsds->getNumberOfDS() == 1, true);

  std::cout << "------- test insertDynamicalSystem ok -------" <<std::endl;
}

// insertInteraction
void NonSmoothDynamicalSystemTest::testinsertInteraction()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(0., 10.));

  SP::DynamicalSystem ds(new LagrangianDS(std::make_shared<SiconosVector>(3),
                                          std::make_shared<SiconosVector>(3)));
  ds->setNumber(23);

  nsds->insertDynamicalSystem(ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertInteractionA: ", nsds->getNumberOfDS() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertInteractionB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertInteractionC: ", nsds->dynamicalSystem(23)->number() == 23, true);

  SP::Relation r(new LagrangianLinearTIR(std::make_shared<SimpleMatrix>(1,3)));
  SP::NonSmoothLaw nsl(new NewtonImpactNSL(0.0));
  SP::Interaction inter(new Interaction(nsl, r));
  nsds->link(inter, ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertInteractionD: ", nsds->getNumberOfInteractions() == 1, true);

  std::cout << "------- test insertInteraction ok -------" <<std::endl;
}


void NonSmoothDynamicalSystemTest::testremoveDynamicalSystem()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(0., 10.));

  SP::DynamicalSystem ds1(new LagrangianDS(std::make_shared<SiconosVector>(3),
                          std::make_shared<SiconosVector>(3)));
  ds1->setNumber(23);
  SP::DynamicalSystem ds2(new LagrangianDS(std::make_shared<SiconosVector>(3),
                          std::make_shared<SiconosVector>(3)));
  ds2->setNumber(32);

  nsds->insertDynamicalSystem(ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemA: ", nsds->getNumberOfDS() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemC: ", nsds->dynamicalSystem(23)->number() == 23, true);

  SP::Relation r1(new LagrangianLinearTIR(std::make_shared<SimpleMatrix>(1,3)));
  SP::Relation r2(new LagrangianLinearTIR(std::make_shared<SimpleMatrix>(1,6)));
  SP::NonSmoothLaw nsl(new NewtonImpactNSL(0.0));
  SP::Interaction inter1(new Interaction(nsl, r1));
  SP::Interaction inter2(new Interaction(nsl, r1));
  SP::Interaction inter3(new Interaction(nsl, r2));
  nsds->link(inter1, ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemD: ", nsds->getNumberOfInteractions() == 1, true);

  nsds->removeDynamicalSystem(ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemE: ", nsds->getNumberOfDS() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemF: ", nsds->getNumberOfInteractions() == 0, true);

  nsds->insertDynamicalSystem(ds1);
  nsds->insertDynamicalSystem(ds2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemG: ", nsds->getNumberOfDS() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemH: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemI: ", nsds->dynamicalSystem(23)->number() == 23, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemI: ", nsds->dynamicalSystem(32)->number() == 32, true);

  nsds->link(inter1, ds1);
  nsds->link(inter2, ds2);
  nsds->link(inter3, ds1, ds2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemJ: ", nsds->getNumberOfInteractions() == 3, true);

  nsds->removeDynamicalSystem(ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemK: ", nsds->getNumberOfDS() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemL: ", nsds->getNumberOfInteractions() == 1, true);

  std::cout << "------- test removeDynamicalSystem ok -------" <<std::endl;
}

void NonSmoothDynamicalSystemTest::testremoveInteraction()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(0., 10.));

  SP::DynamicalSystem ds(new LagrangianDS(std::make_shared<SiconosVector>(3),
                                          std::make_shared<SiconosVector>(3)));
  ds->setNumber(23);

  nsds->insertDynamicalSystem(ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveInteractionA: ", nsds->getNumberOfDS() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveInteractionB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveInteractionC: ", nsds->dynamicalSystem(23)->number() == 23, true);

  SP::Relation r(new LagrangianLinearTIR(std::make_shared<SimpleMatrix>(1,3)));
  SP::NonSmoothLaw nsl(new NewtonImpactNSL(0.0));
  SP::Interaction inter(new Interaction(nsl, r));
  nsds->link(inter, ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveInteractionD: ", nsds->getNumberOfInteractions() == 1, true);

  nsds->removeInteraction(inter);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveInteractionE: ", nsds->getNumberOfInteractions() == 0, true);

  std::cout << "------- test removeInteraction ok -------" <<std::endl;
}
