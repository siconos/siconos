// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "FirstOrderLinearTIRTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearTIRTest);

void FirstOrderLinearTIRTest::setUp() {
  C = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::io::readDenseMatrix("matC.dat"));
  D = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::io::readDenseMatrix("matD.dat"));
  B = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::io::readDenseMatrix("matB.dat"));
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
