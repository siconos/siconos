// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __QuaternionTest__
#define __QuaternionTest__

#include <cppunit/extensions/HelperMacros.h>

#include "RotationQuaternion.hpp"
#include "SiconosVector.hpp"

class QuaternionTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(QuaternionTest);

  // tests to be done ...
  CPPUNIT_TEST(testQuaternion);
  CPPUNIT_TEST(testQuaternionMatrix);
  CPPUNIT_TEST_SUITE_END();

  void testQuaternion();
  void testQuaternionMatrix();

  // Members
 public:
  void setUp();
  void tearDown();
  siconos::algebra::SiconosVector7 q0 = {1, 2, 3, 0, 1, 0, 0};
  siconos::algebra::SiconosVector7 q01 = {1, 2, 3, 1, 0, 0, 0};
};

#endif
