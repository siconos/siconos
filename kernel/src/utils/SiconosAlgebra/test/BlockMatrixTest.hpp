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
#ifndef __BlockMatrixTest__
#define __BlockMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosPointers.hpp"
#include "BlockMatrix.hpp"



class BlockMatrixTest : public CppUnit::TestFixture
{


private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BlockMatrixTest);


  // test suite
  CPPUNIT_TEST_SUITE(BlockMatrixTest);

  CPPUNIT_TEST(testConstructor0);
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testEye);
  CPPUNIT_TEST(testNormInf);
  CPPUNIT_TEST(testGetSetRowCol);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testOperators1);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_SUITE_END();

  void testConstructor0();
  void testConstructor1();
  void testConstructor2();
  void testConstructor3();
  void testConstructor4();
  void testGetSetRowCol();
  void testZero();
  void testEye();
  void testNormInf();
  void testAssignment();
  void testOperators1();
  void End();

  SP::SiconosMatrix B, C, D, E, F, G;
  std::vector<SP::SiconosMatrix> m;
  Index tRow, tCol;
  SP::BlocksMat mapRef;
  double tol;

public:
  void setUp();
  void tearDown();

};

#endif
