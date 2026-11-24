/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __NewtonEulerDSTest__
#define __NewtonEulerDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
class NewtonEulerDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(NewtonEulerDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildNewtonEulerDS1_alias);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildNewtonEulerDS1_alias();
  // void testcomputeDS();

  // Members
  siconos::algebra::SiconosVector q0{7}, q01{7}, twist0{6};
  siconos::algebra::SiconosMatrix inertia{3, 3};
  double mass{10.};

 public:
  void setUp();
  void tearDown();
};

#endif
