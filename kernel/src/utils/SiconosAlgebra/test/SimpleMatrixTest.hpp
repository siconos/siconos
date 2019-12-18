/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef __SimpleMatrixTest__
#define __SimpleMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "BlockMatrix.hpp"
#include "BlockVector.hpp"

class SimpleMatrixTest : public CppUnit::TestFixture
{


private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SimpleMatrixTest);


  // test suite
  CPPUNIT_TEST_SUITE(SimpleMatrixTest);

  CPPUNIT_TEST(testConstructor0);
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testConstructor5);
  CPPUNIT_TEST(testConstructor6);
  CPPUNIT_TEST(testConstructor7);
  CPPUNIT_TEST(testConstructor8);
  CPPUNIT_TEST(testConstructor9);
  CPPUNIT_TEST(testConstructor10);
  CPPUNIT_TEST(testConstructor11);
  CPPUNIT_TEST(testConstructor12);
  CPPUNIT_TEST(testConstructor13);
  CPPUNIT_TEST(testConstructor14);
  CPPUNIT_TEST(testGetSetRowCol);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testEye);
  CPPUNIT_TEST(testResize);
  CPPUNIT_TEST(testNormInf);
  CPPUNIT_TEST(testSetBlock);
  CPPUNIT_TEST(testSetBlock2);
  CPPUNIT_TEST(testTrans);
  CPPUNIT_TEST(testAssignment0);
  CPPUNIT_TEST(testAssignment1);
  CPPUNIT_TEST(testAssignment2);
  CPPUNIT_TEST(testOperators1);
  CPPUNIT_TEST(testOperators2);
  CPPUNIT_TEST(testOperators3);
  CPPUNIT_TEST(testOperators4);
  CPPUNIT_TEST(testOperators4bis);
  CPPUNIT_TEST(testOperators5);
  CPPUNIT_TEST(testOperators6);
  CPPUNIT_TEST(testOperators6Bis);
  CPPUNIT_TEST(testOperators6Ter);
  CPPUNIT_TEST(testOperators7);
  //  CPPUNIT_TEST(testOperators8);
  CPPUNIT_TEST(testOperators8Bis);
  CPPUNIT_TEST(testOperators8Ter);
  CPPUNIT_TEST(testOperators8_4);
  CPPUNIT_TEST(testOperators8_5);
  CPPUNIT_TEST(testOperators8_6);
  CPPUNIT_TEST(testOperators9);
  CPPUNIT_TEST(testOperators9Bis);
  CPPUNIT_TEST(testOperators9Ter);
  CPPUNIT_TEST(testOperators10);
  CPPUNIT_TEST(testOperators11);
  CPPUNIT_TEST(testOperators12);
  CPPUNIT_TEST(testOperators13);
  //CPPUNIT_TEST(testPow);
  CPPUNIT_TEST(testProd);
  CPPUNIT_TEST(testProdBis);
  CPPUNIT_TEST(testProdTer);
  CPPUNIT_TEST(testProd4);
  CPPUNIT_TEST(testProd5);
  CPPUNIT_TEST(testProd6);
  // CPPUNIT_TEST(testGemv);
  // CPPUNIT_TEST(testGemm);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_SUITE_END();

  void testConstructor0();
  void testConstructor1();
  void testConstructor2();
  void testConstructor3();
  void testConstructor4();
  void testConstructor5();
  void testConstructor6();
  void testConstructor7();
  void testConstructor8();
  void testConstructor9();
  void testConstructor10();
  void testConstructor11();
  void testConstructor12();
  void testConstructor13();
  void testConstructor14();
  void testGetSetRowCol();
  void testZero();
  void testEye();
  void testResize();
  void testNormInf();
  void testSetBlock();
  void testSetBlock2();
  void testTrans();
  void testAssignment0();
  void testAssignment1();
  void testAssignment2();
  void testOperators1();
  void testOperators2();
  void testOperators3();
  void testOperators4();
  void testOperators4bis();
  void testOperators5();
  void testOperators6();
  void testOperators6Bis();
  void testOperators6Ter();
  void testOperators7();
  // void testOperators8();
  void testOperators8Bis();
  void testOperators8Ter();
  void testOperators8_4();
  void testOperators8_5();
  void testOperators8_6();
  void testOperators9();
  void testOperators9Bis();
  void testOperators9Ter();
  void testOperators10();
  void testOperators11();
  void testOperators12();
  void testOperators13();
  void testProd();
  void testProdBis();
  void testProdTer();
  void testProd4();
  void testProd5();
  void testProd6();
  // void testGemm();
  // void testGemv();
  void End();

  unsigned int size, size2;
  SP::SiconosMatrix SicM, m1, m2, m3, m4, m5, m6, m7, m8, C, Cb, Cb2;
  SPC::SiconosMatrix A, B, Ab, Bb;
  SP::SimpleMatrix SimM;
  std::string fic1, fic2;
  SP::SiconosVector vect1, vect2, vect3;
  SP::DenseMat  D;
  SP::TriangMat T, T2;
  SP::SymMat S, S2;
  SP::BandedMat Band, Band2;
  SP::SparseMat SP, SP2;
  SP::SparseCoordinateMat SP_coor;
  SP::ZeroMat  Z, Z2;
  SP::IdentityMat I, I2;
  double tol;

public:
  void setUp();
  void tearDown();

};

#endif
