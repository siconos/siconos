/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __RelationCompatibilityTest__
#define __RelationCompatibilityTest__

#include <cppunit/extensions/HelperMacros.h>

/** Checks the deprecated names of the LagrangianR / NewtonEulerR API (kept for one release):
 *  - a class overriding a former virtual name (e.g. computeJacobianhOver_q) is still called
 *    through the new name (computeContactMatrix), which is the one used by the kernel;
 *  - a class overriding the new name is reached through the former name;
 *  - without override, both names run the default implementation (user function);
 *  - the deprecated accessors, setters, members and class names forward to the new ones.
 */
class RelationCompatibilityTest : public CppUnit::TestFixture {
 private:
  CPPUNIT_TEST_SUITE(RelationCompatibilityTest);

  CPPUNIT_TEST(testLagrangianScleronomousR);
  CPPUNIT_TEST(testLagrangianRheonomousR);
  CPPUNIT_TEST(testLagrangianCompliantR);
  CPPUNIT_TEST(testNewtonEulerR);
  CPPUNIT_TEST(testClassAliases);
  CPPUNIT_TEST_SUITE_END();

  void testLagrangianScleronomousR();
  void testLagrangianRheonomousR();
  void testLagrangianCompliantR();
  void testNewtonEulerR();
  void testClassAliases();

 public:
  void setUp() {}
  void tearDown() {}
};

#endif
