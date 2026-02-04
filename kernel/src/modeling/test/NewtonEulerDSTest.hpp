/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef __NewtonEulerDSTest__
#define __NewtonEulerDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
class NewtonEulerDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(NewtonEulerDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildNewtonEulerDS1_alias);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildNewtonEulerDS1_alias();
  // void testcomputeDS();

  // Members
  siconos::algebra::SiconosVector q0{7}, q01{7}, twist0{6};
  siconos::algebra::SiconosMatrix inertia{3, 3};
  double mass{10.};

 public:
  void setUp();
  void tearDown();
};

#endif
