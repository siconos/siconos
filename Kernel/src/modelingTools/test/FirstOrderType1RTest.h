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
#ifndef __FirstOrderType1RTest__
#define __FirstOrderType1RTest__

#include <cppunit/extensions/HelperMacros.h>
#include "RelationXML.h"
#include "NonSmoothDynamicalSystem.h"
#include "NonSmoothDynamicalSystemXML.h"
#include "FirstOrderType1R.h"

class FirstOrderType1RTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderType1RTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderType1R1);
  CPPUNIT_TEST(testBuildFirstOrderType1R2);
  CPPUNIT_TEST(testBuildFirstOrderType1R3);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderType1R1();
  void testBuildFirstOrderType1R2();
  void testBuildFirstOrderType1R3();
  void End();

  // Members

  xmlNodePtr node1;
  SP::RelationXML tmpxml1;
  SP::NonSmoothDynamicalSystem nsds;

public:
  void setUp();
  void tearDown();

};

#endif




