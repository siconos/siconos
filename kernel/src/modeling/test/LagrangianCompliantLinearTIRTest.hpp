// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __LagrangianCompliantLinearTIRTest__
#define __LagrangianCompliantLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianCompliantLinearTIR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class LagrangianCompliantLinearTIRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianCompliantLinearTIRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianCompliantLinearTIRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianCompliantLinearTIR1);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianCompliantLinearTIR1();
  void End();

  // Members

  std::shared_ptr<siconos::algebra::SiconosMatrix> C, B, F, D;
  std::shared_ptr<siconos::algebra::SiconosVector> e;
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds;

 public:
  void setUp();
  void tearDown();
};

#endif
