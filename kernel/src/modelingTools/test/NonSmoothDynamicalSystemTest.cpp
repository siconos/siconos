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

#include "Interaction.hpp"
#include "LagrangianLinearTIR.hpp"
#include "NewtonImpactNSL.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
#include "Topology.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NonSmoothDynamicalSystemTest);

void NonSmoothDynamicalSystemTest::setUp() {}

void NonSmoothDynamicalSystemTest::tearDown() {}

// insertDynamicalSystem
void NonSmoothDynamicalSystemTest::testinsertDynamicalSystem()
{
  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(1., 10.);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" test ndsn build: ", nsds->currentTime() == 1., true);

  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(
      std::make_shared<siconos::algebra::SiconosVector>(3),
      std::make_shared<siconos::algebra::SiconosVector>(3));
  ds->setNumber(23);

  try {
    std::shared_ptr<siconos::modeling::DynamicalSystem> dsnull;
    nsds->insertDynamicalSystem(dsnull);
  }
  catch (const siconos::exception& e) {
    /*  Pass */
    std::cout << "testinsertDynamicalSystemNull: success!" << std::endl;
  }
  catch (const std::exception& e) {
    std::cout << "testinsertDynamicalSystemNull:" << e.what() << std::endl;
    CPPUNIT_FAIL("testinsertDynamicalSystemNull: unexpected exception ");
  }

  nsds->insertDynamicalSystem(ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemA: ", nsds->getNumberOfDS() == 1,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testinsertDynamicalSystemB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testinsertDynamicalSystemC: ", nsds->dynamicalSystem(23)->number() == 23, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testinsertDynamicalSystemD: ", nsds->topology()->hasDynamicalSystem(ds) == true, true);

  // Try again : must be ignored
  nsds->insertDynamicalSystem(ds);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemE: ", nsds->getNumberOfDS() == 1,
                               true);

  std::cout << "------- test insertDynamicalSystem ok -------" << std::endl;
}

// insertInteraction
void NonSmoothDynamicalSystemTest::testinsertInteraction()
{
  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0., 10.);

  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(
      std::make_shared<siconos::algebra::SiconosVector>(3),
      std::make_shared<siconos::algebra::SiconosVector>(3));
  ds->setNumber(23);

  nsds->insertDynamicalSystem(ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertInteractionA: ", nsds->getNumberOfDS() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testinsertInteractionB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testinsertInteractionC: ", nsds->dynamicalSystem(23)->number() == 23, true);

  auto r = std::make_shared<siconos::modeling::LagrangianLinearTIR>(
      std::make_shared<siconos::algebra::SiconosMatrix>(1, 3));
  auto nsl = std::make_shared<siconos::modeling::NewtonImpactNSL>(0.0);
  auto inter = std::make_shared<siconos::modeling::Interaction>(nsl, r);

  nsds->link(inter, ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testinsertInteractionD: ", nsds->getNumberOfInteractions() == 1, true);

  std::cout << "------- test insertInteraction ok -------" << std::endl;
}

void NonSmoothDynamicalSystemTest::testremoveDynamicalSystem()
{
  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0., 10.);

  auto ds1 = std::make_shared<siconos::modeling::LagrangianDS>(
      std::make_shared<siconos::algebra::SiconosVector>(3),
      std::make_shared<siconos::algebra::SiconosVector>(3));

  ds1->setNumber(23);
  auto ds2 = std::make_shared<siconos::modeling::LagrangianDS>(
      std::make_shared<siconos::algebra::SiconosVector>(3),
      std::make_shared<siconos::algebra::SiconosVector>(3));
  ds2->setNumber(32);

  nsds->insertDynamicalSystem(ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemA: ", nsds->getNumberOfDS() == 1,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemC: ", nsds->dynamicalSystem(23)->number() == 23, true);

  auto r1 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(
      std::make_shared<siconos::algebra::SiconosMatrix>(1, 3));
  auto r2 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(
      std::make_shared<siconos::algebra::SiconosMatrix>(1, 6));
  auto nsl = std::make_shared<siconos::modeling::NewtonImpactNSL>(0.0);
  auto inter1 = std::make_shared<siconos::modeling::Interaction>(nsl, r1);
  auto inter2 = std::make_shared<siconos::modeling::Interaction>(nsl, r1);
  auto inter3 = std::make_shared<siconos::modeling::Interaction>(nsl, r2);
  nsds->link(inter1, ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemD: ", nsds->getNumberOfInteractions() == 1, true);

  nsds->removeDynamicalSystem(ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemE: ", nsds->getNumberOfDS() == 0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemF: ", nsds->getNumberOfInteractions() == 0, true);

  nsds->insertDynamicalSystem(ds1);
  nsds->insertDynamicalSystem(ds2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemG: ", nsds->getNumberOfDS() == 2,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemH: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemI: ", nsds->dynamicalSystem(23)->number() == 23, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemI: ", nsds->dynamicalSystem(32)->number() == 32, true);

  nsds->link(inter1, ds1);
  nsds->link(inter2, ds2);
  nsds->link(inter3, ds1, ds2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemJ: ", nsds->getNumberOfInteractions() == 3, true);

  nsds->removeDynamicalSystem(ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveDynamicalSystemK: ", nsds->getNumberOfDS() == 1,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveDynamicalSystemL: ", nsds->getNumberOfInteractions() == 1, true);

  std::cout << "------- test removeDynamicalSystem ok -------" << std::endl;
}

void NonSmoothDynamicalSystemTest::testremoveInteraction()
{
  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0., 10.);

  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(
      std::make_shared<siconos::algebra::SiconosVector>(3),
      std::make_shared<siconos::algebra::SiconosVector>(3));
  ds->setNumber(23);

  nsds->insertDynamicalSystem(ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testremoveInteractionA: ", nsds->getNumberOfDS() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveInteractionB: ", nsds->getNumberOfInteractions() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveInteractionC: ", nsds->dynamicalSystem(23)->number() == 23, true);

  auto r = std::make_shared<siconos::modeling::LagrangianLinearTIR>(
      std::make_shared<siconos::algebra::SiconosMatrix>(1, 3));
  auto nsl = std::make_shared<siconos::modeling::NewtonImpactNSL>(0.0);
  auto inter = std::make_shared<siconos::modeling::Interaction>(nsl, r);

  nsds->link(inter, ds);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveInteractionD: ", nsds->getNumberOfInteractions() == 1, true);

  nsds->removeInteraction(inter);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      " testremoveInteractionE: ", nsds->getNumberOfInteractions() == 0, true);

  std::cout << "------- test removeInteraction ok -------" << std::endl;
}
