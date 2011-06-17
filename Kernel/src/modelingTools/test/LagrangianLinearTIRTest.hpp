/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
#ifndef __LagrangianLinearTIRTest__
#define __LagrangianLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "RelationXML.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NonSmoothDynamicalSystemXML.hpp"
#include "LagrangianLinearTIR.hpp"
#include "LinearRXML.hpp"

class LagrangianLinearTIRTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianLinearTIRTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianLinearTIRTest);

  // tests to be done ...

  //CPPUNIT_TEST(testBuildLagrangianLinearTIR);
  CPPUNIT_TEST(testBuildLagrangianLinearTIR0);
  CPPUNIT_TEST(testBuildLagrangianLinearTIR1);
  CPPUNIT_TEST(testBuildLagrangianLinearTIR2);
  CPPUNIT_TEST(testBuildLagrangianLinearTIR3);
  CPPUNIT_TEST(testBuildLagrangianLinearTIR4);
  CPPUNIT_TEST(testBuildLagrangianLinearTIR5);
  CPPUNIT_TEST(testBuildLagrangianLinearTIR6);
  CPPUNIT_TEST(testSetCPtr);
  CPPUNIT_TEST(testSetDPtr);
  CPPUNIT_TEST(testSetFPtr);
  CPPUNIT_TEST(testSetEPtr);
  CPPUNIT_TEST(testGetJacPtr);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianLinearTIR0();
  void testBuildLagrangianLinearTIR1();
  void testBuildLagrangianLinearTIR2();
  void testBuildLagrangianLinearTIR3();
  void testBuildLagrangianLinearTIR4();
  void testBuildLagrangianLinearTIR5();
  void testBuildLagrangianLinearTIR6();
  void testSetCPtr();
  void testSetDPtr();
  void testSetFPtr();
  void testSetEPtr();
  void testGetJacPtr();
  void End();

  // Members

  SP::SiconosMatrix C, B, F, D;
  SP::SiconosVector e;
  xmlNodePtr node1;
  SP::RelationXML tmpxml1;
  SP::NonSmoothDynamicalSystem nsds;

public:
  void setUp();
  void tearDown();

};

#endif




