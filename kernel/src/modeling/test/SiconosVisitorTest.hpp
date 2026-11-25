// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef SiconosVisitorTest_h
#define SiconosVisitorTest_h

#include <cppunit/extensions/HelperMacros.h>

class SiconosVisitorTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SiconosVisitorTest);

  // tests to be done ...
  CPPUNIT_TEST(testVisitNSL);
  CPPUNIT_TEST_SUITE_END();

  // Members
  void testVisitNSL();

 public:
  void setUp();
  void tearDown();
};

#endif
