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
#ifndef __CABLEDSTest__
#define __CABLEDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include "TransportCableManager.h"

class CableDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(CableDSTest);

  CPPUNIT_TEST(testReadModel);
  CPPUNIT_TEST(testBuildInitialProfile);
  CPPUNIT_TEST(testComputeDS);
  CPPUNIT_TEST(testComputeBouncingBall);
  CPPUNIT_TEST(testNoFext);
  CPPUNIT_TEST(testConstantFext);
  CPPUNIT_TEST(testVariableFext);

  CPPUNIT_TEST_SUITE_END();

  void testReadModel();
  void testBuildInitialProfile();
  void testComputeDS();
  void testComputeBouncingBall();

  void testNoFext();
  void testConstantFext();
  void testVariableFext();

 public:
  void setUp();
  void tearDown();
};

#endif
