// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __LagrangianScleronomousRTest__
#define __LagrangianScleronomousRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianScleronomousR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class LagrangianScleronomousRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianScleronomousRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianScleronomousRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianScleronomousR2);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianScleronomousR0();
  void testBuildLagrangianScleronomousR2();

 public:
  void setUp();
  void tearDown();
};

#endif
