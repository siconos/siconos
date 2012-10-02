/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#ifndef __FirstOrderType1RTest__
#define __FirstOrderType1RTest__

#include <cppunit/extensions/HelperMacros.h>
#include "RelationXML.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NonSmoothDynamicalSystemXML.hpp"
#include "FirstOrderType1R.hpp"

class FirstOrderType1RTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderType1RTest);


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




