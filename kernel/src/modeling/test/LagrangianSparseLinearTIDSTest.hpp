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
#ifndef __LagrangianSparseLinearTIDSTest__
#define __LagrangianSparseLinearTIDSTest__

#include <cppunit/extensions/HelperMacros.h>

#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>

class LagrangianSparseLinearTIDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianSparseLinearTIDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianSparseLinearTIDS_basic);
  CPPUNIT_TEST(testBuildLagrangianSparseLinearTIDS_alias);
  CPPUNIT_TEST(testBuildLagrangianSparseLinearTIDS_copy);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianSparseLinearTIDS_basic();
  void testBuildLagrangianSparseLinearTIDS_alias();
  void testBuildLagrangianSparseLinearTIDS_copy();
  // void testcomputeDS();

  // Members

  siconos::algebra::SiconosVector3 q0, velocity0;

  std::shared_ptr<siconos::algebra::SiconosSparseMatrix> mass{nullptr};

 public:
  void setUp();
  void tearDown();
};

#endif
