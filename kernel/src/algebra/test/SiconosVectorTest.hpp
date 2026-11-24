// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>

#include "SiconosVector.hpp"

class SiconosVectorTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SiconosVectorTest);

  CPPUNIT_TEST(testOrthoBaseFromVector);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();
  void testOrthoBaseFromVector();
  void End();

 public:
  void setUp();
  void tearDown();
};

#endif
