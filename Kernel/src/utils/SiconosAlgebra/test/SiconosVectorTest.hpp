/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>
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




