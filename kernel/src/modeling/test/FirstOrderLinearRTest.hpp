// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __FirstOrderLinearRTest__
#define __FirstOrderLinearRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "FirstOrderLinearR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class FirstOrderLinearRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(FirstOrderLinearRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearR1);
  CPPUNIT_TEST(testBuildFirstOrderLinearR2);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearR1();
  void testBuildFirstOrderLinearR2();
  // Members

  std::shared_ptr<siconos::algebra::SiconosMatrix> C, B, D;
  std::shared_ptr<siconos::algebra::SiconosVector> e;
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds;

 public:
  void setUp();
  void tearDown();
};

#endif
