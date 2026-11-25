// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __FirstOrderNonLinearDSTest__
#define __FirstOrderNonLinearDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class FirstOrderNonLinearDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderNonLinearDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS_alias);
  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS_copy);
  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS2);
  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS3);
  CPPUNIT_TEST(testInitMemory);
  CPPUNIT_TEST(testSwap);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderNonLinearDS_alias();
  void testBuildFirstOrderNonLinearDS_copy();
  void testBuildFirstOrderNonLinearDS2();
  void testBuildFirstOrderNonLinearDS3();
  void testInitMemory();
  void testSwap();

  // Members

  std::shared_ptr<siconos::algebra::SiconosVector> x0{nullptr}, xnull{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> J0{nullptr}, M{nullptr};

 public:
  void setUp();
  void tearDown();
};

#endif
