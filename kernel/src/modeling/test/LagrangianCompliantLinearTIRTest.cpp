// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LagrangianCompliantLinearTIRTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianCompliantLinearTIRTest);

void LagrangianCompliantLinearTIRTest::setUp() {
  C = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::io::readDenseMatrix("matC.dat"));
  D = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::io::readDenseMatrix("matD.dat"));
  e = std::make_shared<siconos::algebra::SiconosVector>(1);
  (*e)(0) = 0.1;
}

void LagrangianCompliantLinearTIRTest::tearDown() {}

// data constructor (1)
void LagrangianCompliantLinearTIRTest::testBuildLagrangianCompliantLinearTIR1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto folr = std::make_shared<siconos::modeling::LagrangianCompliantLinearTIR>(*C, *D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianCompliantLinearTIR1a : ", folr->contactMatrix() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR1a : ", folr->D() == *D,
                               true);
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
      "testBuildLagrangianCompliantLinearTIR1a : ", folr2->contactMatrix() == *C, true);
  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

void LagrangianCompliantLinearTIRTest::End() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== End of LagrangianCompliantLinearTIR Tests ===== " << std::endl;
  std::cout << "=========================================== " << std::endl;
}
