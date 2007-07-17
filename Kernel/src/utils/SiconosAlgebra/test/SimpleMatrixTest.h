/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#ifndef __SimpleMatrixTest__
#define __SimpleMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SimpleMatrix.h"
#include "BlockMatrix.h"
#include "BlockVector.h"

using namespace std;

class SimpleMatrixTest : public CppUnit::TestFixture
{


private:

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
  CPPUNIT_TEST(testOperators5);
  CPPUNIT_TEST(testOperators6);
  CPPUNIT_TEST(testOperators6Bis);
  CPPUNIT_TEST(testOperators6Ter);
  CPPUNIT_TEST(testOperators7);
  CPPUNIT_TEST(testOperators8);
  CPPUNIT_TEST(testOperators8Bis);
  CPPUNIT_TEST(testOperators8Ter);
  CPPUNIT_TEST(testOperators8_4);
  CPPUNIT_TEST(testOperators9);
  CPPUNIT_TEST(testOperators9Bis);
  CPPUNIT_TEST(testOperators9Ter);
  CPPUNIT_TEST(testOperators10);
  CPPUNIT_TEST(testOperators11);
  CPPUNIT_TEST(testOperators12);
  CPPUNIT_TEST(testOperators13);
  CPPUNIT_TEST(testPow);
  CPPUNIT_TEST(testProd);
  CPPUNIT_TEST(testProdBis);
  CPPUNIT_TEST(testProdTer);
  CPPUNIT_TEST(testProd4);
  CPPUNIT_TEST(testProd5);
  CPPUNIT_TEST(testProd6);
  CPPUNIT_TEST(testGemv);
  CPPUNIT_TEST(testGemm);
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
  void testOperators5();
  void testOperators6();
  void testOperators6Bis();
  void testOperators6Ter();
  void testOperators7();
  void testOperators8();
  void testOperators8Bis();
  void testOperators8Ter();
  void testOperators8_4();
  void testOperators9();
  void testOperators9Bis();
  void testOperators9Ter();
  void testOperators10();
  void testOperators11();
  void testOperators12();
  void testOperators13();
  void testPow();
  void testProd();
  void testProdBis();
  void testProdTer();
  void testProd4();
  void testProd5();
  void testProd6();
  void testGemm();
  void testGemv();
  void End();

  unsigned int size, size2;
  SiconosMatrix *SicM, *m1, *m2, *m3, *m4, *m5, *m6, *m7, *m8, *C, *Cb, *Cb2;
  const SiconosMatrix *A, *B, *Ab, *Bb;
  SimpleMatrix *SimM;
  string fic1, fic2;
  SimpleVector* vect1, *vect2, *vect3;
  DenseMat * D;
  TriangMat *T, *T2;
  SymMat *S, *S2;
  BandedMat *Band, *Band2;
  SparseMat *SP;
  ZeroMat * Z, *Z2;
  IdentityMat* I, *I2;
  double tol;

public:
  void setUp();
  void tearDown();

};

#endif
