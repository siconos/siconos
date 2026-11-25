// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __CABLEDSTest__
#define __CABLEDSTest__

#include <cppunit/extensions/HelperMacros.h>

class CableDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(CableDSTest);

  CPPUNIT_TEST(testReadModel);
  CPPUNIT_TEST(testBuildInitialProfile);
  CPPUNIT_TEST(testInitializeFEM);
  CPPUNIT_TEST(testComputeDS);
  CPPUNIT_TEST(testComputeBouncingBall);
  CPPUNIT_TEST(testNoFext);
  CPPUNIT_TEST(testConstantFext);
  CPPUNIT_TEST(testVariableFext);

  CPPUNIT_TEST_SUITE_END();

  void testReadModel();
  void testBuildInitialProfile();
  void testInitializeFEM();
  void testComputeDS();
  void testComputeBouncingBall();

  void testNoFext();
  void testConstantFext();
  void testVariableFext();

 public:
  void setUp();
  void tearDown();
};

#endif
