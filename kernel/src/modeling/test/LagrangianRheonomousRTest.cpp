// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LagrangianRheonomousRTest.hpp"

#include "Interaction.hpp"
#include "NewtonImpactNSL.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianRheonomousRTest);

void LagrangianRheonomousRTest::setUp() {}

void LagrangianRheonomousRTest::tearDown() {}

// Complete test of LagrangianRheonomousR: see CamFollower example.

//  constructor:
void LagrangianRheonomousRTest::testBuildLagrangianRheonomousR0() {
  auto rel = std::make_shared<siconos::modeling::LagrangianRheonomousR>();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianRheonomousR3a : ",
                               rel->getType() == siconos::modeling::RelationType::Lagrangian,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianRheonomousR3b : ",
      rel->getSubType() == siconos::modeling::RelationSubType::RheonomousR, true);

  auto hfunc = [](const siconos::algebra::BlockVector& pos, double time,
                  Eigen::Ref<siconos::algebra::MapVectorType> result) { result.setZero(); };

  auto jachq = [](const siconos::algebra::BlockVector& pos, double time,
                  Eigen::Ref<siconos::algebra::MapType> result) { result.setZero(); };

  auto dhdt = [](const siconos::algebra::BlockVector& pos, double time,
                 Eigen::Ref<siconos::algebra::MapVectorType> result) { result.setZero(); };

  rel->setComputehFunction(hfunc);

  rel->setComputeContactMatrixFunction(jachq);

  rel->setComputedhdtFunction(dhdt);

  // For a complete test/example with initialize(inter), compute call and so on
  // see CamFollower-Rheonomous
  std::cout << " data Constructor LagrangianRheonomousR ok" << std::endl;
}
