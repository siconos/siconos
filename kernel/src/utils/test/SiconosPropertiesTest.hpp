// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef SiconosPropertiesTest_h
#define SiconosPropertiesTest_h

#include <cppunit/extensions/HelperMacros.h>

#include "SiconosGraph.hpp"
#include "SiconosProperties.hpp"

class SiconosPropertiesTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SiconosPropertiesTest);

  // tests to be done ...
  CPPUNIT_TEST(t1);

  CPPUNIT_TEST(t2);

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
