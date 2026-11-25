/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __FirstOrderLinearTIRTest__
#define __FirstOrderLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "FirstOrderLinearTIR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class FirstOrderLinearTIRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(FirstOrderLinearTIRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearTIRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearTIR1);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearTIR1();

  // Members

  std::shared_ptr<siconos::algebra::SiconosMatrix> C, B, D;
  std::shared_ptr<siconos::algebra::SiconosVector> e;

 public:
  void setUp();
  void tearDown();
};

#endif
