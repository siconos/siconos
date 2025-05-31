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
#ifndef __FirstOrderLinearTIRTest__
#define __FirstOrderLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>

#include "FirstOrderLinearTIR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class FirstOrderLinearTIRTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(FirstOrderLinearTIRTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearTIRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearTIR1);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearTIR1();

  // Members

  std::shared_ptr<siconos::algebra::SiconosMatrix> C, B, D;
  std::shared_ptr<siconos::algebra::SiconosVector> e;

 public:
  void setUp();
  void tearDown();
};

#endif
