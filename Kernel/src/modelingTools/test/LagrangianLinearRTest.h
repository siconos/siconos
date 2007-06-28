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
#ifndef __LagrangianLinearRTest__
#define __LagrangianLinearRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LagrangianLinearR.h"
#include "LagrangianLinearRXML.h"
#include "NonSmoothDynamicalSystem.h"
#include "NonSmoothDynamicalSystemXML.h"

class LagrangianLinearRTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianLinearRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianLinearR0);
  CPPUNIT_TEST(testBuildLagrangianLinearR1);
  CPPUNIT_TEST(testBuildLagrangianLinearR2);
  CPPUNIT_TEST(testBuildLagrangianLinearR1);
  CPPUNIT_TEST(testBuildLagrangianLinearR4);
  CPPUNIT_TEST(testSetH);
  CPPUNIT_TEST(testSetHPtr);
  CPPUNIT_TEST(testSetB);
  CPPUNIT_TEST(testSetBPtr);
  CPPUNIT_TEST(testSetD);
  CPPUNIT_TEST(testSetDPtr);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianLinearR0();
  void testBuildLagrangianLinearR1();
  void testBuildLagrangianLinearR2();
  void testBuildLagrangianLinearR3();
  void testBuildLagrangianLinearR4();
  void testSetH();
  void testSetHPtr();
  void testSetB();
  void testSetBPtr();
  void testSetD();
  void testSetDPtr();
  void End();

  // Members

  SiconosMatrix *H, *D, *F;
  SimpleVector *b;
  LagrangianRXML* tmpxml1;

public:
  void setUp();
  void tearDown();

};

#endif




