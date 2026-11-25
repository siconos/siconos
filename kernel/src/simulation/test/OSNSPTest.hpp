// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __OSNSPTest__
#define __OSNSPTest__

#include <cppunit/extensions/HelperMacros.h>

class OSNSPTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(OSNSPTest);

  // tests to be done ...
  CPPUNIT_TEST(testOSNSBuild_default);
  CPPUNIT_TEST(testOSNSBuild_solverid);
  CPPUNIT_TEST(testOSNSBuild_options);
  CPPUNIT_TEST_SUITE_END();

  void testOSNSBuild_default();
  void testOSNSBuild_solverid();
  void testOSNSBuild_options();

 public:
  void setUp();
  void tearDown();
};

#endif
