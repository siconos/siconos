// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LagrangianLinearTIRTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearTIRTest);

void LagrangianLinearTIRTest::setUp() {
  C = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::io::readDenseMatrix("matC.dat"));
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
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIR1a : ", folr->contactMatrix() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1c : ",
                               folr->getType() == siconos::modeling::RelationType::Lagrangian,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIR1d : ",
      folr->getSubType() == siconos::modeling::RelationSubType::LinearTIR, true);

  auto folr2 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*C, *e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1a : ", folr2->eVector() == *e,
                               true);

  (*e)(0) = 4.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1a : ", folr2->eVector() == *e,
                               true);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}
