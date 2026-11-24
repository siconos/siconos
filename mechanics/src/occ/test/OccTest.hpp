/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef OccTest_h
#define OccTest_h

#include <cppunit/extensions/HelperMacros.h>

class OccTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(OccTest);

  CPPUNIT_TEST(exportBRepAsString);

  CPPUNIT_TEST(computeUVBounds);

  CPPUNIT_TEST(move);
#ifdef HAS_FORTRAN
  CPPUNIT_TEST(distance);
#endif
  CPPUNIT_TEST_SUITE_END();

  // Members
  void exportBRepAsString();

  void computeUVBounds();

  void move();

#ifdef HAS_FORTRAN
  void distance();
#endif

 public:
  void setUp();
  void tearDown();
};

#endif
