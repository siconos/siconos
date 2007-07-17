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
#ifndef __FirstOrderNonLinearDSTest__
#define __FirstOrderNonLinearDSTest__

#include <cppunit/extensions/HelperMacros.h>
#include "FirstOrderNonLinearDSXML.h"
#include "FirstOrderNonLinearDS.h"
#include "RuntimeException.h"
#include "XMLException.h"

class FirstOrderNonLinearDSTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderNonLinearDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS1);
  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS2);
  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS3);
  CPPUNIT_TEST(testSetX0);
  CPPUNIT_TEST(testSetX0Ptr);
  CPPUNIT_TEST(testSetX);
  CPPUNIT_TEST(testSetXPtr);
  CPPUNIT_TEST(testSetXFree);
  CPPUNIT_TEST(testSetXFreePtr);
  CPPUNIT_TEST(testSetR);
  CPPUNIT_TEST(testSetRPtr);
  CPPUNIT_TEST(testSetJacobianXF);
  CPPUNIT_TEST(testSetJacobianXFPtr);
  CPPUNIT_TEST(testInitMemory);
  CPPUNIT_TEST(testSwap);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderNonLinearDS1();
  void testBuildFirstOrderNonLinearDS2();
  void testBuildFirstOrderNonLinearDS3();
  void testSetX0();
  void testSetX0Ptr();
  void testSetX();
  void testSetXPtr();
  void testSetXFree();
  void testSetXFreePtr();
  void testSetR();
  void testSetRPtr();
  void testSetJacobianXF();
  void testSetJacobianXFPtr();
  void testInitMemory();
  void testSwap();
  void End();

  // Members

  SimpleVector * x0;
  SiconosMatrix *J0, *M;
  xmlNode * node1 , *node2;
  FirstOrderNonLinearDSXML* tmpxml1, * tmpxml2;
public:
  void setUp();
  void tearDown();

};

#endif




