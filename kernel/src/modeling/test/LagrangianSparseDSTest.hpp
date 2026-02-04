/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#ifndef __LagrangianSparseDSTest__
#define __LagrangianSparseDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include "LagrangianSparseDS.hpp"
#include "SiconosException.hpp"

class LagrangianSparseDSTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(LagrangianSparseDSTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianSparseDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianSparseDS_basic);
  CPPUNIT_TEST(testBuildLagrangianSparseDS_alias);
  CPPUNIT_TEST(testBuildLagrangianSparseDS_copy);
  CPPUNIT_TEST(testBuildLagrangianSparseDS_compute);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianSparseDS_basic();
  void testBuildLagrangianSparseDS_alias();
  void testBuildLagrangianSparseDS_copy();
  void testBuildLagrangianSparseDS_compute();
  // void testcomputeDS();

  // Members

  siconos::algebra::SiconosVector3 q0, velocity0;

  std::shared_ptr<siconos::algebra::SiconosSparseMatrix> mass{nullptr};

 public:
  void setUp();
  void tearDown();
};

#endif
