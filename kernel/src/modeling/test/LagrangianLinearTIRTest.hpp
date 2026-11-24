// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __LagrangianLinearTIRTest__
#define __LagrangianLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianLinearTIR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class LagrangianLinearTIRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianLinearTIRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianLinearTIRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianLinearTIR1);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianLinearTIR1();

  // Members

  std::shared_ptr<siconos::algebra::SiconosMatrix> C;
  std::shared_ptr<siconos::algebra::SiconosVector> e;

 public:
  void setUp();
  void tearDown();
};

#endif
