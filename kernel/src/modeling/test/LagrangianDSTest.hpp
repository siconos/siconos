// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __LagrangianDSTest__
#define __LagrangianDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianDS.hpp"
#include "SiconosException.hpp"

class LagrangianDSTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianDSTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianDS_copy);
  CPPUNIT_TEST(testBuildLagrangianDS_alias);
  CPPUNIT_TEST(testBuildLagrangianDS2_alias);
  CPPUNIT_TEST(testBuildLagrangianDS2_copy);
  CPPUNIT_TEST(testBuildLagrangianDS3);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianDS_copy();
  void testBuildLagrangianDS_alias();
  void testBuildLagrangianDS2_alias();
  void testBuildLagrangianDS2_copy();
  void testBuildLagrangianDS3();
  // void testcomputeDS();

  // Members

  siconos::algebra::SiconosVector3 q0, velocity0;

  siconos::algebra::SiconosMatrix mass{3, 3};

 public:
  void setUp();
  void tearDown();
};

#endif
