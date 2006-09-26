/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
#ifndef __MySimpleMatrixTest__
#define __MySimpleMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "MySimpleMatrix.h"

using namespace std;

class MySimpleMatrixTest : public CppUnit::TestFixture
{


private:

  // test suite
  CPPUNIT_TEST_SUITE(MySimpleMatrixTest);

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
  CPPUNIT_TEST(testGetSetRowCol);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testEye);
  CPPUNIT_TEST(testResize);
  CPPUNIT_TEST(testNormInf);
  CPPUNIT_TEST(testGetBlock);
  CPPUNIT_TEST(testBlockMatrixCopy);
  CPPUNIT_TEST(testOperators);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testMultTranspose);
  CPPUNIT_TEST(End);

  /*   CPPUNIT_TEST(testReadWriteBinary); */
  /*   CPPUNIT_TEST(testReadWriteAscii); */
  /*   CPPUNIT_TEST(testAssignment); */
  /*   CPPUNIT_TEST(testMultTranspose); */
  /*   CPPUNIT_TEST(testBlockMatrixCopy1); */
  /*   CPPUNIT_TEST(testBlockMatrixCopy2); */
  /*   CPPUNIT_TEST(testBlockMatrixCopy3); */
  /*   CPPUNIT_TEST(testBlockMatrixCopy4); */
  /*   CPPUNIT_TEST(End); */
  /*   CPPUNIT_TEST_EXCEPTION(testBlockMatrixCopyException1, SiconosMatrixException); */
  /*   CPPUNIT_TEST_EXCEPTION(testBlockMatrixCopyException2, SiconosMatrixException); */
  /*   CPPUNIT_TEST(testOperators); */
  /*   CPPUNIT_TEST(testLinearSolve); */
  /*   CPPUNIT_TEST_EXCEPTION(testSizeException, SiconosMatrixException); */
  /*   CPPUNIT_TEST_EXCEPTION(testConstructorException, SiconosMatrixException); */

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
  void testGetSetRowCol();
  void testZero();
  void testEye();
  void testResize();
  void testNormInf();
  void testGetBlock();
  void testBlockMatrixCopy();
  void testOperators();
  void testAssignment();
  void testMultTranspose();
  void End();


  /*   void testReadWriteBinary(); */
  /*   void testReadWriteAscii(); */
  /*   void testBlockMatrixCopy1(); */
  /*   void testBlockMatrixCopy2(); */
  /*   void testBlockMatrixCopy3(); */
  /*   void testBlockMatrixCopy4(); */
  /*   void testBlockMatrixCopyException1(); */
  /*   void testBlockMatrixCopyException2(); */
  /*   void testLinearSolve(); */
  /*   void testSizeException(); */
  /*   void testConstructorException(); */

  MySiconosMatrix *SicM;
  MySimpleMatrix *SimM;
  string fic1, fic2;
  MySimpleVector* vect1, *vect2, *vect3;

public:
  void setUp();
  void tearDown();

};

#endif
