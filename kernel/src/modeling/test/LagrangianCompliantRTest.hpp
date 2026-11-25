/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __LagrangianCompliantRTest__
#define __LagrangianCompliantRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianCompliantR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class LagrangianCompliantRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianCompliantRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianCompliantRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianCompliantR0);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianCompliantR0();
  // Members

 public:
  void setUp();
  void tearDown();
};

#endif
