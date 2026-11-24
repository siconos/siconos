/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __LagrangianRheonomousRTest__
#define __LagrangianRheonomousRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianRheonomousR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class LagrangianRheonomousRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianRheonomousRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianRheonomousRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianRheonomousR0);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianRheonomousR0();

 public:
  void setUp();
  void tearDown();
};

#endif
