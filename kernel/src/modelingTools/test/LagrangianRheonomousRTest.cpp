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

  auto hfunc = [](const siconos::algebra::BlockVector &pos, double time,
                  Eigen::Ref<siconos::algebra::MapVectorType> result) { result.setZero(); };

  auto jachq = [](const siconos::algebra::BlockVector &pos, double time,
                  Eigen::Ref<siconos::algebra::MapType> result) { result.setZero(); };

  auto hdot = [](const siconos::algebra::BlockVector &pos, double time,
                 Eigen::Ref<siconos::algebra::MapVectorType> result) { result.setZero(); };

  rel->setComputehFunction(hfunc);

  rel->setComputeJacobianhOver_qFunction(jachq);

  rel->setComputehdotFunction(hdot);

  // For a complete test/example with initialize(inter), compute call and so on
  // see CamFollower-Rheonomous
  std::cout << " data Constructor LagrangianRheonomousR ok" << std::endl;
}
