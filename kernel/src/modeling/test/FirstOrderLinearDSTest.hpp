/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __FirstOrderLinearDSTest__
#define __FirstOrderLinearDSTest__
#include <cppunit/extensions/HelperMacros.h>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class FirstOrderLinearDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearDS1_alias);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS1_copy);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS2_alias);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS2_copy);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS_plugins);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearDS1_alias();
  void testBuildFirstOrderLinearDS1_copy();
  void testBuildFirstOrderLinearDS2_alias();
  void testBuildFirstOrderLinearDS2_copy();
  void testBuildFirstOrderLinearDS_plugins();

  // Members

  std::shared_ptr<siconos::algebra::SiconosVector> x0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> b0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> A0{nullptr};

 public:
  void setUp();
  void tearDown();
};

#endif
