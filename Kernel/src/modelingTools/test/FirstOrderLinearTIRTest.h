/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#ifndef __FirstOrderLinearTIRTest__
#define __FirstOrderLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "RelationXML.h"
#include "NonSmoothDynamicalSystem.h"
#include "NonSmoothDynamicalSystemXML.h"
#include "FirstOrderLinearTIR.h"
#include "FirstOrderLinearTIRXML.h"

class FirstOrderLinearTIRTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearTIRTest);

  // tests to be done ...

  //CPPUNIT_TEST(testBuildFirstOrderLinearTIR);
  CPPUNIT_TEST(testBuildFirstOrderLinearTIR0);
  CPPUNIT_TEST(testBuildFirstOrderLinearTIR1);
  CPPUNIT_TEST(testBuildFirstOrderLinearTIR2);
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
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearTIR0();
  void testBuildFirstOrderLinearTIR1();
  void testBuildFirstOrderLinearTIR2();
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
  void End();

  // Members

  SiconosMatrix *C, *B, *F, *D;
  SimpleVector *e;
  xmlNode * node1;
  RelationXML * tmpxml1;
  NonSmoothDynamicalSystem * nsds;

public:
  void setUp();
  void tearDown();

};

#endif




