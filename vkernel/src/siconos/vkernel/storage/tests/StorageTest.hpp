/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <cppunit/extensions/HelperMacros.h>

class StorageTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(StorageTest);
  CPPUNIT_TEST(testStorage0);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp() override;
  void tearDown() override;
  void testStorage0();
};
