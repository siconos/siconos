/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
#ifndef __LinearTIRTest__
#define __LinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LinearTIR.h"

class LinearTIRTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LinearTIRTest);

  // tests to be done ...

  //CPPUNIT_TEST(testBuildLinearTIR);
  CPPUNIT_TEST(testBuildLinearTIR0);
  CPPUNIT_TEST(testBuildLinearTIR1);
  CPPUNIT_TEST(testBuildLinearTIR2);
  CPPUNIT_TEST(testBuildLinearTIR3);
  CPPUNIT_TEST(testBuildLinearTIR4);
  CPPUNIT_TEST(testSetC);
  CPPUNIT_TEST(testSetCPtr);
  CPPUNIT_TEST(testSetD);
  CPPUNIT_TEST(testSetDPtr);
  CPPUNIT_TEST(testSetF);
  CPPUNIT_TEST(testSetFPtr);
  CPPUNIT_TEST(testSetE);
  CPPUNIT_TEST(testSetEPtr);
  CPPUNIT_TEST(testSetB);
  CPPUNIT_TEST(testSetBPtr);
  CPPUNIT_TEST(testSetA);
  CPPUNIT_TEST(testSetAPtr);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLinearTIR0();
  void testBuildLinearTIR1();
  void testBuildLinearTIR2();
  void testBuildLinearTIR3();
  void testBuildLinearTIR4();
  void testSetC();
  void testSetCPtr();
  void testSetD();
  void testSetDPtr();
  void testSetF();
  void testSetFPtr();
  void testSetE();
  void testSetEPtr();
  void testSetB();
  void testSetBPtr();
  void testSetA();
  void testSetAPtr();
  void End();

  // Members

  SiconosMatrix *C, *B, *F, *D;
  SimpleVector *a, *e;
  xmlNode * node1, *node2;
  RelationXML * tmpxml1, *tmpxml2;
  NonSmoothDynamicalSystem * nsds;

public:
  void setUp();
  void tearDown();

};

#endif




