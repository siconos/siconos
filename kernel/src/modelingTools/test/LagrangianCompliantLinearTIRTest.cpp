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
#include "LagrangianCompliantLinearTIRTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianCompliantLinearTIRTest);

void LagrangianCompliantLinearTIRTest::setUp() {
  C = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matC.dat"));
  D = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matD.dat"));
  e = std::make_shared<siconos::algebra::SiconosVector>(1);
  (*e)(0) = 0.1;
}

void LagrangianCompliantLinearTIRTest::tearDown() {}

// data constructor (1)
void LagrangianCompliantLinearTIRTest::testBuildLagrangianCompliantLinearTIR1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto folr = std::make_shared<siconos::modeling::LagrangianCompliantLinearTIR>(*C, *D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantLinearTIR1a : ", folr->CMatrix() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantLinearTIR1a : ", folr->DMatrix() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR1c : ",
                               folr->getType() == siconos::modeling::RelationType::Lagrangian,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantLinearTIR1d : ",
      folr->getSubType() == siconos::modeling::RelationSubType::CompliantLinearTIR, true);

  auto folr2 = std::make_shared<siconos::modeling::LagrangianCompliantLinearTIR>(*C, *D, *e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantLinearTIR1a : ", folr2->eVector() == *e, true);

  (*e)(0) = 4.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantLinearTIR1a : ", folr2->eVector() == *e, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantLinearTIR1a : ", folr2->jacobianhOver_q() == *C, true);
  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

void LagrangianCompliantLinearTIRTest::End() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== End of LagrangianCompliantLinearTIR Tests ===== " << std::endl;
  std::cout << "=========================================== " << std::endl;
}
