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
#ifndef __FirstOrderLinearRTest__
#define __FirstOrderLinearRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "FirstOrderLinearR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class FirstOrderLinearRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(FirstOrderLinearRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearR1);
  CPPUNIT_TEST(testBuildFirstOrderLinearR2);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearR1();
  void testBuildFirstOrderLinearR2();
  // Members

  std::shared_ptr<siconos::algebra::SiconosMatrix> C, B, D;
  std::shared_ptr<siconos::algebra::SiconosVector> e;
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds;

 public:
  void setUp();
  void tearDown();
};

#endif
