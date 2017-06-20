/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#ifndef __BlockVectorTest__
#define __BlockVectorTest__

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <cppunit/extensions/HelperMacros.h>
#include "SiconosVector.hpp"
#include <cmath>
#include <vector>

class BlockVectorTest : public CppUnit::TestFixture
{


private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BlockVectorTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(BlockVectorTest);

  // tests to be done ...

  //  CPPUNIT_TEST(testBuildBlockVector);
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testConstructor5);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testFill);
  CPPUNIT_TEST(testNorm);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testOperators1);
  CPPUNIT_TEST(testOperators2);
  CPPUNIT_TEST(testOperators3);
  CPPUNIT_TEST(testOperators4);
  CPPUNIT_TEST(testInsert);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testConstructor1();
  void testConstructor2();
  void testConstructor3();
  void testConstructor4();
  void testConstructor5();
  void testZero();
  void testFill();
  void testNorm();
  void testAssignment();
  void testOperators1();
  void testOperators2();
  void testOperators3();
  void testOperators4();
  void testInsert();
  void End();
  // Members

  SP::BlockVector ref;
  std::vector<double> vq;
  SP::DenseVect dv;
  SP::SparseVect sv;
  double tol;

public:
  void setUp();
  void tearDown();

};

#endif




