// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LagrangianScleronomousRTest.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianScleronomousRTest);

void LagrangianScleronomousRTest::setUp() {}

void LagrangianScleronomousRTest::tearDown() {}

// data constructor:
void LagrangianScleronomousRTest::testBuildLagrangianScleronomousR2() {
  auto rel = std::make_shared<siconos::modeling::LagrangianScleronomousR>();

  // Just set user-defined functions operators.
  // We can not call compute... functions since in relations
  // everything is properly set (memory) only after a call to initialize
  // which required an Interaction. See examples in siconos tutorials for reals tests.

  rel->setComputehFunction([](const siconos::algebra::BlockVector& q,
                              Eigen::Ref<siconos::algebra::SiconosVector> y) {});

  rel->setComputeContactMatrixFunction([](const siconos::algebra::BlockVector& q,
                                          Eigen::Ref<siconos::algebra::MapType> result) {});

  rel->setComputeContactMatrixDotFunction([](const siconos::algebra::BlockVector& q,
                                             const siconos::algebra::BlockVector& qdot,
                                             Eigen::Ref<siconos::algebra::MapType> result) {});

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR3a : ",
                               rel->getType() == siconos::modeling::RelationType::Lagrangian,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianScleronomousR3b : ",
      rel->getSubType() == siconos::modeling::RelationSubType::ScleronomousR, true);
  std::cout << " data Constructor LagrangianScleronomousR ok" << std::endl;
}
