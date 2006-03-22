/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#ifndef __BlockMatrixTest__
#define __BlockMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "BlockMatrix.h"
#include "SimpleMatrix.h"

using namespace std;

class BlockMatrixTest : public CppUnit::TestFixture
{

private:

  // test suite
  CPPUNIT_TEST_SUITE(BlockMatrixTest);

  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testConstructor3);
  CPPUNIT_TEST(testConstructor4);
  CPPUNIT_TEST(testGetLaGenMatDouble);
  CPPUNIT_TEST(testSetValue);
  CPPUNIT_TEST(testAssignment);
  CPPUNIT_TEST(testOperators);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  void testConstructor1();
  void testConstructor2();
  void testConstructor3();
  void testConstructor4();
  void testGetLaGenMatDouble();
  void testSetValue();
  void testAssignment();
  void testOperators();
  void End();

  SiconosMatrix *B, *C, *D, *E, *F, *G;
  vector<SiconosMatrix*> m;

public:
  void setUp();
  void tearDown();

};

#endif
