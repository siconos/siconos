/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __LagrangianSparseDSTest__
#define __LagrangianSparseDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class LagrangianSparseDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianSparseDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianSparseDS_basic);
  CPPUNIT_TEST(testBuildLagrangianSparseDS_alias);
  CPPUNIT_TEST(testBuildLagrangianSparseDS_copy);
  CPPUNIT_TEST(testBuildLagrangianSparseDS_compute);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianSparseDS_basic();
  void testBuildLagrangianSparseDS_alias();
  void testBuildLagrangianSparseDS_copy();
  void testBuildLagrangianSparseDS_compute();
  // void testcomputeDS();

  // Members

  siconos::algebra::SiconosVector3 q0, velocity0;

  std::shared_ptr<siconos::algebra::SiconosSparseMatrix> mass{nullptr};

 public:
  void setUp();
  void tearDown();
};

#endif
