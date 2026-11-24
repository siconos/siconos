/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __SiconosMatrixTest__
#define __SiconosMatrixTest__

#include <cppunit/extensions/HelperMacros.h>

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class SiconosMatrixTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(SiconosMatrixTest);

  // test suite
  CPPUNIT_TEST_SUITE(SiconosMatrixTest);

  CPPUNIT_TEST(testNormInf);
  CPPUNIT_TEST(testSymm);
  CPPUNIT_TEST(testfillTriplet);
  CPPUNIT_TEST(testBuildSparse);
  CPPUNIT_TEST(testSetBlock);
  CPPUNIT_TEST(testProd);
  CPPUNIT_TEST(testProd2);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_SUITE_END();

  void testNormInf();
  void testSymm();
  void testfillTriplet();
  void testBuildSparse();
  void testSetBlock();
  void testProd();
  void testProd2();

  void End();

  unsigned int size;
  std::shared_ptr<siconos::algebra::SiconosDenseMatrix> SicM;
  siconos::algebra::SiconosSparseMatrix Asparse;
  std::string fic1, fic2;
  double tol;

 public:
  void setUp();
  void tearDown();
};

#endif
