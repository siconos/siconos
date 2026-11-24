/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __NonSmoothDynamicalSystemTest__
#define __NonSmoothDynamicalSystemTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianLinearTIDS.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class NonSmoothDynamicalSystemTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(NonSmoothDynamicalSystemTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(NonSmoothDynamicalSystemTest);

  // tests to be done ...

  // CPPUNIT_TEST(testBuildNonSmoothDynamicalSystem);
  CPPUNIT_TEST(testinsertDynamicalSystem);
  CPPUNIT_TEST(testinsertInteraction);
  CPPUNIT_TEST(testremoveDynamicalSystem);
  CPPUNIT_TEST(testremoveInteraction);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testinsertDynamicalSystem();
  void testinsertInteraction();
  void testremoveDynamicalSystem();
  void testremoveInteraction();

 public:
  void setUp();
  void tearDown();
};

#endif
