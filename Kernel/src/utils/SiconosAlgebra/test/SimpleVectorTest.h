/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#ifndef __SimpleVectorTest__
#define __SimpleVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosVector.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include <math.h>
#include <vector>

class SimpleVectorTest : public CppUnit::TestFixture
{


private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SimpleVectorTest);

  // tests to be done ...

  //  CPPUNIT_TEST(testBuildSimpleVector);
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testConstructor5);
  CPPUNIT_TEST(testConstructor6);
  CPPUNIT_TEST(testConstructor7);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testNorm);
  CPPUNIT_TEST(testResize);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testOperators1);
  CPPUNIT_TEST(testOperators2);
  CPPUNIT_TEST(testOperators3);
  CPPUNIT_TEST(testOperators4);
  CPPUNIT_TEST(testOperators5);
  CPPUNIT_TEST(testOperators6);
  CPPUNIT_TEST(testOperators7);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testConstructor1();
  void testConstructor2();
  void testConstructor3();
  void testConstructor4();
  void testConstructor5();
  void testConstructor6();
  void testConstructor7();
  void testZero();
  void testNorm();
  void testResize();
  void testAssignment();
  void testOperators1();
  void testOperators2();
  void testOperators3();
  void testOperators4();
  void testOperators5();
  void testOperators6();
  void testOperators7();
  void End();
  // Members

  SiconosVector * ref;
  std::vector<double> vq;
  DenseVect * dv;
  SparseVect * sv;
  double tol;

public:
  void setUp();
  void tearDown();

};

#endif




