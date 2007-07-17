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
#ifndef __BlockVectorTest__
#define __BlockVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosVector.h"
#include <math.h>
#include <vector>

class BlockVectorTest : public CppUnit::TestFixture
{


private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(BlockVectorTest);

  // tests to be done ...

  //  CPPUNIT_TEST(testBuildBlockVector);
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
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




