/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#ifndef __SiconosMemoryTest__
#define __SiconosMemoryTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosVector.h"
#include "SiconosMemory.h"
#include "SiconosMemoryException.h"
#include <deque>


#define SIZE 10
#define M_SIZE 3

using namespace std;

class SiconosMemoryTest : public CppUnit::TestFixture
{


private:

  // Test suite
  CPPUNIT_TEST_SUITE(SiconosMemoryTest);

  CPPUNIT_TEST(testBuildMemory1);
  CPPUNIT_TEST(testBuildMemory2);
  CPPUNIT_TEST(testBuildMemory3);
  CPPUNIT_TEST(testBuildMemory4);
  CPPUNIT_TEST(testSetVectorMemory);
  CPPUNIT_TEST(testGetSiconosVector);
  CPPUNIT_TEST(testSwap);
  CPPUNIT_TEST(testOperatorEqual);

  CPPUNIT_TEST_EXCEPTION(testMemoryException, SiconosMemoryException);
  CPPUNIT_TEST_EXCEPTION(testMemoryException1, SiconosMemoryException);

  //CPPUNIT_TEST_FAIL(testFail);

  CPPUNIT_TEST_SUITE_END();

  void testBuildMemory1();
  void testBuildMemory2();
  void testBuildMemory3();
  void testBuildMemory4();
  void testSetVectorMemory();
  void testGetSiconosVector();
  void testSwap();
  void testOperatorEqual();

  void testMemoryException();
  void testMemoryException1();
  //void testFail();

  deque<SiconosVector*> V1;
  deque<SiconosVector*> V2;
  deque<SiconosVector*> V3;
  SiconosVector * q1;
  SiconosVector * q2;
  SiconosVector * q3;
  SiconosVector *c1;
  SiconosVector *c2;
  unsigned int sizeMem;
public:
  void setUp();
  void tearDown();

};

#endif

