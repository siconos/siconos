// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __LagrangianSparseLinearTIDSTest__
#define __LagrangianSparseLinearTIDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>

class LagrangianSparseLinearTIDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianSparseLinearTIDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianSparseLinearTIDS_basic);
  CPPUNIT_TEST(testBuildLagrangianSparseLinearTIDS_alias);
  CPPUNIT_TEST(testBuildLagrangianSparseLinearTIDS_copy);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianSparseLinearTIDS_basic();
  void testBuildLagrangianSparseLinearTIDS_alias();
  void testBuildLagrangianSparseLinearTIDS_copy();
  // void testcomputeDS();

  // Members

  siconos::algebra::SiconosVector3 q0, velocity0;

  std::shared_ptr<siconos::algebra::SiconosSparseMatrix> mass{nullptr};

 public:
  void setUp();
  void tearDown();
};

#endif
