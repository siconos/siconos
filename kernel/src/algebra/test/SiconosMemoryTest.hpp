// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __SiconosMemoryTest__
#define __SiconosMemoryTest__

#include <cppunit/extensions/HelperMacros.h>

#include "BlockVector.hpp"
#include "SiconosMemory.hpp"

class SiconosMemoryTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(SiconosMemoryTest);

  // Test suite
  CPPUNIT_TEST_SUITE(SiconosMemoryTest);

  CPPUNIT_TEST(testBuildMemory1);
  // CPPUNIT_TEST(testBuildMemory2);
  // CPPUNIT_TEST(testBuildMemory3);
  CPPUNIT_TEST(testSetVectorMemory);
  CPPUNIT_TEST(testGetSiconosVector);
  CPPUNIT_TEST(testSwap);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_SUITE_END();

  void testBuildMemory1();
  // void testBuildMemory2();
  // void testBuildMemory3();
  void testSetVectorMemory();
  void testGetSiconosVector();
  void testSwap();
  void End();

  std::shared_ptr<siconos::algebra::blocks::Vector> V1, V2, V3;
  std::shared_ptr<siconos::algebra::SiconosVector> q1, q2, q3;
  std::shared_ptr<siconos::algebra::BlockVector> c1, c2;
  unsigned int _sizeMem;

 public:
  void setUp();
  void tearDown();
};

#endif
