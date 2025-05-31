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
#include "FirstOrderLinearTIRTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearTIRTest);

void FirstOrderLinearTIRTest::setUp() {
  C = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matC.dat"));
  D = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matD.dat"));
  B = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matB.dat"));
  e = std::make_shared<siconos::algebra::SiconosVector>(1);
  (*e)(0) = 0.1;
}

void FirstOrderLinearTIRTest::tearDown() {}

// data constructor (1)
void FirstOrderLinearTIRTest::testBuildFirstOrderLinearTIR1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto folr = std::make_shared<siconos::modeling::FirstOrderLinearTIR>();
  folr->setConstantB(*B);
  folr->setConstantC(*C);
  folr->setConstantD(*D);
  folr->setConstanteVector(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1a : ", folr->B() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1a : ", folr->C() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1a : ", folr->D() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1a : ", folr->eVector() == *e,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1c : ",
                               folr->getType() == siconos::modeling::RelationType::FirstOrder,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderLinearTIR1d : ",
      folr->getSubType() == siconos::modeling::RelationSubType::LinearTIR, true);
  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}
