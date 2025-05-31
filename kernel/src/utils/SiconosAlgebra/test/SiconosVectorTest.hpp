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
#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <cmath>
#include <vector>

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class SiconosVectorTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(SiconosVectorTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SiconosVectorTest);

  // tests to be done ...

  //  CPPUNIT_TEST(testBuildSiconosVector);
  CPPUNIT_TEST(testSetBlockFriend);
  CPPUNIT_TEST(testSetBlock);
  CPPUNIT_TEST(testOperators4Ter);
  CPPUNIT_TEST(testOperators5Bis);
  CPPUNIT_TEST(testOperators6Bis);
  CPPUNIT_TEST(testOperators7);
  CPPUNIT_TEST(testOrthoBaseFromVector);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testSetBlock();
  void testSetBlockFriend();
  void testOperators4Ter();
  void testOperators5Bis();
  void testOperators6Bis();
  void testOperators7();
  void testOrthoBaseFromVector();
  void End();
  // Members

  std::shared_ptr<siconos::algebra::SiconosVector> ref, z, tmp1, tmp2, tmp3, tmp4;
  std::shared_ptr<siconos::algebra::BlockVector> zB;
  std::shared_ptr<const siconos::algebra::SiconosVector> x, y;
  std::shared_ptr<const siconos::algebra::BlockVector> xB, yB;
  unsigned int size, size1, size2;
  std::vector<double> vq;
  double tol;

 public:
  void setUp();
  void tearDown();
};

#endif
