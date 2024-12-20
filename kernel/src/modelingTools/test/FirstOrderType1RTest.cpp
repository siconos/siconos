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

#include "FirstOrderType1RTest.hpp"

#include "Interaction.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderType1RTest);

void FirstOrderType1RTest::setUp() {}

void FirstOrderType1RTest::tearDown() {}

// data constructor
void FirstOrderType1RTest::testBuildFirstOrderType1R1() {
  std::cout << "======================================" << std::endl;
  std::cout << "=== FirstOrderType1R tests start ...== " << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "--> Test: constructor 1 " << std::endl;

  // Just set constant or user-defined functions operators.
  // We can not call compute... functions since in relations
  // everything is properly set (memory) only after a call to initialize
  // which required an Interaction. See examples in siconos tutorials for reals tests.

  auto rel_cst = std::make_shared<siconos::modeling::FirstOrderType1R>();

  rel_cst->setComputehFunction([](const siconos::algebra::BlockVector &state,
                                  Eigen::Ref<siconos::algebra::SiconosVector> y) {});

  rel_cst->setComputegFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda,
         siconos::algebra::BlockVector &res) {});

  siconos::algebra::SiconosMatrix cst_mat{3, 3};
  cst_mat.setRandom();
  rel_cst->setConstantJacobianhOver_state(cst_mat);

  rel_cst->setConstantJacobiangOver_lambda(cst_mat);

  auto rel = std::make_shared<siconos::modeling::FirstOrderType1R>();

  rel->setComputehFunction([](const siconos::algebra::BlockVector &state,
                              Eigen::Ref<siconos::algebra::SiconosVector> y) {});

  rel->setComputegFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda,
                              siconos::algebra::BlockVector &res) {});

  rel->setComputeJacobianhOver_stateFunction(
      [](const siconos::algebra::BlockVector &state,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  rel->setComputeJacobiangOver_lambdaFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1b : ",
                               rel->getType() == siconos::modeling::RelationType::FirstOrder,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1relc : ",
                               rel->getSubType() == siconos::modeling::RelationSubType::Type1R,
                               true);
  std::cout << "--> Constructor1 test ended with success." << std::endl;
}
