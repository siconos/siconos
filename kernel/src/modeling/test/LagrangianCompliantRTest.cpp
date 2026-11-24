/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "LagrangianCompliantRTest.hpp"

#include "SiconosMatrix.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianCompliantRTest);

void LagrangianCompliantRTest::setUp() {}

void LagrangianCompliantRTest::tearDown() {}

// data constructor:
void LagrangianCompliantRTest::testBuildLagrangianCompliantR0() {
  auto rel = std::make_shared<siconos::modeling::LagrangianCompliantR>();

  // Just set user-defined functions operators.
  // We can not call compute... functions since in relations
  // everything is properly set (memory) only after a call to initialize
  // which required an Interaction. See examples in siconos tutorials for reals tests.

  rel->setComputehFunction([](const siconos::algebra::BlockVector& q,
                              const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                              Eigen::Ref<siconos::algebra::SiconosVector> y) {});

  rel->setComputeContactMatrixFunction(
      [](const siconos::algebra::BlockVector& q,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  rel->setComputeDFunction([](const siconos::algebra::BlockVector& q,
                              const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                              Eigen::Ref<siconos::algebra::MapType> result) {});

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3a : ",
                               rel->getType() == siconos::modeling::RelationType::Lagrangian,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantR3b : ",
      rel->getSubType() == siconos::modeling::RelationSubType::CompliantR, true);
  std::cout << " data Constructor LagrangianCompliantR ok" << std::endl;
}
