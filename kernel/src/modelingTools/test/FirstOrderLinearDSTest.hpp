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
#ifndef __FirstOrderLinearDSTest__
#define __FirstOrderLinearDSTest__
#include <cppunit/extensions/HelperMacros.h>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class FirstOrderLinearDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearDS1_alias);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS1_copy);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS2_alias);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS2_copy);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS_plugins);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearDS1_alias();
  void testBuildFirstOrderLinearDS1_copy();
  void testBuildFirstOrderLinearDS2_alias();
  void testBuildFirstOrderLinearDS2_copy();
  void testBuildFirstOrderLinearDS_plugins();

  // Members

  std::shared_ptr<siconos::algebra::SiconosVector> x0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> b0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> A0{nullptr};

 public:
  void setUp();
  void tearDown();
};

#endif
