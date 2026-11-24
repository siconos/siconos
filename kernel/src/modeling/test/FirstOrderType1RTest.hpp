/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __FirstOrderType1RTest__
#define __FirstOrderType1RTest__

#include <cppunit/extensions/HelperMacros.h>

#include "FirstOrderType1R.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class FirstOrderType1RTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(FirstOrderType1RTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderType1RTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderType1R1);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderType1R1();

  // Members

  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds;

 public:
  void setUp();
  void tearDown();
};

#endif
