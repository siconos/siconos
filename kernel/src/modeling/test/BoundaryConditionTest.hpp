/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BoundaryConditionTest__
#define __BoundaryConditionTest__

#include <cppunit/extensions/HelperMacros.h>

#include "BoundaryCondition.hpp"

class BoundaryConditionTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(BoundaryConditionTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(BoundaryConditionTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildBoundaryConditionFixed);
  CPPUNIT_TEST(testBuildBoundaryConditionWithValues);
  CPPUNIT_TEST(testBuildHarmonicBoundaryConditionScalar);
  CPPUNIT_TEST(testBuildHarmonicBoundaryConditionVector);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildBoundaryConditionFixed();
  void testBuildBoundaryConditionWithValues();
  void testBuildHarmonicBoundaryConditionScalar();
  void testBuildHarmonicBoundaryConditionVector();

 public:
  void setUp() {};
  void tearDown() {};
};

#endif
