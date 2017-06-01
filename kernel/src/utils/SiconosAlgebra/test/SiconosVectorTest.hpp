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
#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include <cmath>
#include <vector>

class SiconosVectorTest : public CppUnit::TestFixture
{


private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosVectorTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SiconosVectorTest);

  // tests to be done ...

  //  CPPUNIT_TEST(testBuildSiconosVector);
  CPPUNIT_TEST(testConstructor0);
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor3Bis);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testConstructor5);
  CPPUNIT_TEST(testConstructor6);
  CPPUNIT_TEST(testConstructor7);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testFill);
  CPPUNIT_TEST(testNorm);
  CPPUNIT_TEST(testResize);
  CPPUNIT_TEST(testSetBlock);
  //  CPPUNIT_TEST(testSetBlock2);
  CPPUNIT_TEST(testSetBlock3);
  CPPUNIT_TEST(testSetBlock4);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testOperators1);
  CPPUNIT_TEST(testOperators2);
  CPPUNIT_TEST(testOperators3);
  CPPUNIT_TEST(testOperators4);
  CPPUNIT_TEST(testOperators4Bis);
  CPPUNIT_TEST(testOperators4Ter);
  CPPUNIT_TEST(testOperators5);
  CPPUNIT_TEST(testOperators5Bis);
  CPPUNIT_TEST(testOperators6);
  CPPUNIT_TEST(testOperators6Bis);
  CPPUNIT_TEST(testOperators7);
  CPPUNIT_TEST(testOperators8);
  CPPUNIT_TEST(testSubscal);
  CPPUNIT_TEST(testStdOstream);
  CPPUNIT_TEST(testStdVectorCast);
  CPPUNIT_TEST(testIterators);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testConstructor0();
  void testConstructor1();
  void testConstructor2();
  void testConstructor3();
  void testConstructor3Bis();
  void testConstructor4();
  void testConstructor5();
  void testConstructor6();
  void testConstructor7();
  void testZero();
  void testFill();
  void testNorm();
  void testResize();
  void testSetBlock();
  void testSetBlock2();
  void testSetBlock3();
  void testSetBlock4();
  void testAssignment();
  void testOperators1();
  void testOperators2();
  void testOperators3();
  void testOperators4();
  void testOperators4Bis();
  void testOperators4Ter();
  void testOperators5();
  void testOperators5Bis();
  void testOperators6();
  void testOperators6Bis();
  void testOperators7();
  void testOperators8();
  void testSubscal();
  void testStdOstream();
  void testStdVectorCast();
  void testIterators();
  void End();
  // Members

  SP::SiconosVector ref, z, tmp1, tmp2, tmp3, tmp4;
  SP::BlockVector zB;
  SPC::SiconosVector x, y;
  SPC::BlockVector xB, yB;
  unsigned int size, size1, size2;
  std::vector<double> vq;
  SP::DenseVect  dv;
  SP::SparseVect  sv;
  double tol;

public:
  void setUp();
  void tearDown();

};

#endif




