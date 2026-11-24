// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MechanicsIO_TEST_HPP
#define MechanicsIO_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>

class MechanicsIOTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(MechanicsIOTest);

  CPPUNIT_TEST(test1);
  CPPUNIT_TEST_SUITE_END();

  void test1();

 public:
  void setUp();
  void tearDown();
};

#endif
