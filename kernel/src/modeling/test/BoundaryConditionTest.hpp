/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef __BoundaryConditionTest__
#define __BoundaryConditionTest__

#include <cppunit/extensions/HelperMacros.h>

#include "BoundaryCondition.hpp"

class BoundaryConditionTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(BoundaryConditionTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(BoundaryConditionTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildBoundaryConditionFixed);
  CPPUNIT_TEST(testBuildBoundaryConditionWithValues);
  CPPUNIT_TEST(testBuildHarmonicBoundaryConditionScalar);
  CPPUNIT_TEST(testBuildHarmonicBoundaryConditionVector);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildBoundaryConditionFixed();
  void testBuildBoundaryConditionWithValues();
  void testBuildHarmonicBoundaryConditionScalar();
  void testBuildHarmonicBoundaryConditionVector();

 public:
  void setUp() {};
  void tearDown() {};
};

#endif
