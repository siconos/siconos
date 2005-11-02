/* Siconos version 1.0, Copyright INRIA 2005.
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
#ifndef __LagrangianRTest__
#define __LagrangianRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LagrangianR.h"

class LagrangianRTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianR0);
  CPPUNIT_TEST(testBuildLagrangianR1);
  CPPUNIT_TEST(testBuildLagrangianR4);
  CPPUNIT_TEST(testBuildLagrangianR2);
  CPPUNIT_TEST(testBuildLagrangianR3);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianR0();
  void testBuildLagrangianR1();
  void testBuildLagrangianR4();
  void testBuildLagrangianR2();
  void testBuildLagrangianR3();
  void End();

  // Members

  SiconosMatrix *G0, *G1;
  NonSmoothDynamicalSystem * nsds ;
  xmlNode * node;
  LagrangianRXML* tmpxml1, * tmpxml2, *tmpxml3;
  Interaction * interaction;

public:
  void setUp();
  void tearDown();

};

#endif




