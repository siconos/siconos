/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __LagrangianSparseRelationsTest__
#define __LagrangianSparseRelationsTest__

#include <cppunit/extensions/HelperMacros.h>

class LagrangianSparseRelationsTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianSparseRelationsTest);

  // tests to be done ...
  CPPUNIT_TEST(testBuildLagrangianSparseRheonomousR);
  CPPUNIT_TEST(testBuildLagrangianSparseScleronomousR);
  CPPUNIT_TEST_SUITE_END();

  void testBuildLagrangianSparseRheonomousR();
  void testBuildLagrangianSparseScleronomousR();

 public:
  void setUp();
  void tearDown();
};

#endif
