/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __LagrangianLinearTIDSTest__
#define __LagrangianLinearTIDSTest__

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianLinearTIDS.hpp"
#include "SiconosException.hpp"

class LagrangianLinearTIDSTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianLinearTIDSTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianLinearTIDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianLinearTIDS_alias);
  CPPUNIT_TEST(testBuildLagrangianLinearTIDS_copy);
  CPPUNIT_TEST(testBuildLagrangianLinearTIDS2);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianLinearTIDS_alias();
  void testBuildLagrangianLinearTIDS_copy();
  void testBuildLagrangianLinearTIDS2();

  // Members

  siconos::algebra::SiconosVector3 q0, velocity0;
  siconos::algebra::SiconosMatrix mass{3, 3}, minus_inv_M{3, 3};
  siconos::algebra::SiconosMatrix K{3, 3}, C{3, 3}, rhsK{3, 3}, rhsC{3, 3};

 public:
  void setUp();
  void tearDown();
};

#endif
