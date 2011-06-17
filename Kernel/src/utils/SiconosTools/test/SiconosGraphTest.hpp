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
#ifndef SiconosGraphTest_h
#define SiconosGraphTest_h

#include <cppunit/extensions/HelperMacros.h>
#include "../SiconosGraph.hpp"

class SiconosGraphTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosGraphTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SiconosGraphTest);

  // tests to be done ...
  CPPUNIT_TEST(t1);

  CPPUNIT_TEST(t2);

  CPPUNIT_TEST(t3);

  CPPUNIT_TEST(t4);

  CPPUNIT_TEST(t5);

  CPPUNIT_TEST(t6);

  CPPUNIT_TEST(t7);
  CPPUNIT_TEST(t8);

  CPPUNIT_TEST_SUITE_END();

  // Members
  void t1();
  void t2();
  void t3();
  void t4();
  void t5();
  void t6();
  void t7();
  void t8();

public:
  void setUp();
  void tearDown();

};

#endif
