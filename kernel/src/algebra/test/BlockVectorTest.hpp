/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BlockVectorTest__
#define __BlockVectorTest__

#include <cppunit/extensions/HelperMacros.h>

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <cmath>
#include <vector>

#include "BlockVector.hpp"
#include "SiconosVector.hpp"

class BlockVectorTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(BlockVectorTest);

  // tests to be done ...

  CPPUNIT_TEST(testDefaultConstructor);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testNorm);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testOperators);
  CPPUNIT_TEST(testToSiconosVector);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testDefaultConstructor();
  void testConstructor2();
  void testConstructor3();
  void testNorm();
  void testAssignment();
  void testOperators();
  void testToSiconosVector();
  void End();
  // Members

  double tol{1.e-14};

 public:
  void setUp();
  void tearDown();
};

#endif
