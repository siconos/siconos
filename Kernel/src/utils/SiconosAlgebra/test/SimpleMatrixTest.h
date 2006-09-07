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
#ifndef __SimpleMatrixTest__
#define __SimpleMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SimpleMatrix.h"

using namespace std;

class SimpleMatrixTest : public CppUnit::TestFixture
{


private:

  // test suite
  CPPUNIT_TEST_SUITE(SimpleMatrixTest);

  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testConstructor5);
  CPPUNIT_TEST(testGetRow);
  CPPUNIT_TEST(testReadWriteBinary);
  CPPUNIT_TEST(testReadWriteAscii);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testMultTranspose);
  CPPUNIT_TEST(testBlockMatrixCopy1);
  CPPUNIT_TEST(testBlockMatrixCopy2);
  CPPUNIT_TEST(testBlockMatrixCopy3);
  CPPUNIT_TEST(testBlockMatrixCopy4);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_EXCEPTION(testBlockMatrixCopyException1, SiconosMatrixException);
  CPPUNIT_TEST_EXCEPTION(testBlockMatrixCopyException2, SiconosMatrixException);
  CPPUNIT_TEST(testOperators);
  CPPUNIT_TEST(testLinearSolve);
  CPPUNIT_TEST_EXCEPTION(testSizeException, SiconosMatrixException);
  CPPUNIT_TEST_EXCEPTION(testConstructorException, SiconosMatrixException);

  CPPUNIT_TEST_SUITE_END();

  void testConstructor1();
  void testConstructor2();
  void testConstructor3();
  void testConstructor4();
  void testConstructor5();
  void testGetRow();
  void testReadWriteBinary();
  void testReadWriteAscii();
  void testAssignment();
  void testMultTranspose();
  void testBlockMatrixCopy1();
  void testBlockMatrixCopy2();
  void testBlockMatrixCopy3();
  void testBlockMatrixCopy4();
  void testBlockMatrixCopyException1();
  void testBlockMatrixCopyException2();
  void testOperators();
  void testLinearSolve();
  void testSizeException();
  void testConstructorException();
  void End();

  SiconosMatrix *A, *B, *C;
  SimpleVector* SV;
  LaVectorDouble LVD;

public:
  void setUp();
  void tearDown();

};

#endif
