// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef Test2dContact_h
#define Test2dContact_h

#include <cppunit/extensions/HelperMacros.h>

class Contact2dTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(Contact2dTest);

  // tests to be done ...
  CPPUNIT_TEST(t1);
  // tests to be done ...
  CPPUNIT_TEST(t2);
  // tests to be done ...
  CPPUNIT_TEST(t3);

  CPPUNIT_TEST_SUITE_END();

  // Members
  void t1();
  void t2();
  void t3();

 public:
  void setUp();
  void tearDown();
};

#endif
