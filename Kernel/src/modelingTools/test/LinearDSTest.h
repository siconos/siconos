/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#ifndef __LinearDSTest__
#define __LinearDSTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LinearDS.h"
#include "RuntimeException.h"
#include "XMLException.h"

class LinearDSTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LinearDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLinearDS1);
  CPPUNIT_TEST(testBuildLinearDS2);
  CPPUNIT_TEST(testBuildLinearDS3);
  CPPUNIT_TEST(testBuildLinearDS4);
  CPPUNIT_TEST(testSetA);
  CPPUNIT_TEST(testSetAPtr);
  CPPUNIT_TEST(testSetB);
  CPPUNIT_TEST(testSetBPtr);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLinearDS1();
  void testBuildLinearDS2();
  void testBuildLinearDS3();
  void testBuildLinearDS4();
  void testSetA();
  void testSetAPtr();
  void testSetB();
  void testSetBPtr();
  void End();

  // Members

  SimpleVector * x0, *b0;
  SiconosMatrix *A0, *J0;
  xmlNode * node1 , *node2;
  LinearDSXML* tmpxml1, * tmpxml2;
public:
  void setUp();
  void tearDown();

};

#endif




