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
#ifndef __MyBlockMatrixTest__
#define __MyBlockMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "MyBlockMatrix.h"

using namespace std;

class MyBlockMatrixTest : public CppUnit::TestFixture
{


private:

  // test suite
  CPPUNIT_TEST_SUITE(MyBlockMatrixTest);

  CPPUNIT_TEST(testConstructor0);
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testEye);
  CPPUNIT_TEST(testNormInf);
  CPPUNIT_TEST(testGetSetRowCol);
  CPPUNIT_TEST(testGetBlock);
  CPPUNIT_TEST(testMatrixCopy);
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
  void testGetBlock();
  void testMatrixCopy();
  void testAssignment();
  void testOperators1();
  void End();

  MySiconosMatrix *B, *C, *D, *E, *F, *G;
  std::vector<MySiconosMatrix*> m;
  std::vector<unsigned int> tRow, tCol;
  mapped * mapRef;

public:
  void setUp();
  void tearDown();

};

#endif
