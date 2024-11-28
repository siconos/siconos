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
#include "LagrangianLinearTIRTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearTIRTest);

void LagrangianLinearTIRTest::setUp() {
  C = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matC.dat"));
  e = std::make_shared<siconos::algebra::SiconosVector>(1);
  (*e)(0) = 0.1;
}

void LagrangianLinearTIRTest::tearDown() {}

// data constructor (1)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR1() {
  std::cout << "--> Test: constructor 1." << std::endl;

  // Just set constant or user-defined functions operators.
  // We can not call compute... functions since in relations
  // everything is properly set (memory) only after a call to initialize
  // which required an Interaction. See examples in siconos tutorials for reals tests.

  auto folr = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1a : ", folr->CMatrix() == *C,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1c : ",
                               folr->getType() == siconos::modeling::RelationType::Lagrangian,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIR1d : ",
      folr->getSubType() == siconos::modeling::RelationSubType::LinearTIR, true);

  folr->seteVector(*e);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}
