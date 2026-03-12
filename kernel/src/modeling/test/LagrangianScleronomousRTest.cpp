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

  rel->setComputehFunction([](const siconos::algebra::BlockVector &q,
                              Eigen::Ref<siconos::algebra::SiconosVector> y) {});

  rel->setComputeJacobianhOver_qFunction([](const siconos::algebra::BlockVector &q,
                                            Eigen::Ref<siconos::algebra::MapType> result) {});

  rel->setComputejacobianhOver_q_dotFunction(
      [](const siconos::algebra::BlockVector &q, const siconos::algebra::BlockVector &qdot,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR3a : ",
                               rel->getType() == siconos::modeling::RelationType::Lagrangian,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianScleronomousR3b : ",
      rel->getSubType() == siconos::modeling::RelationSubType::ScleronomousR, true);
  std::cout << " data Constructor LagrangianScleronomousR ok" << std::endl;
}
