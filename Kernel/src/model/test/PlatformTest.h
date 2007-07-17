/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
//$id$

#ifndef PLATFORMTEST_H
#define PLATFORMTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "Model.h"
#include "LCP.h"
#include "RuntimeException.h"
#include "XMLException.h"

using namespace std;

class PlatformTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PlatformTest);

  CPPUNIT_TEST(testManualCreation1);
  CPPUNIT_TEST(testManualCreation2);
  CPPUNIT_TEST(testManualCreation3);
  CPPUNIT_TEST_EXCEPTION(testPlatformException, RuntimeException);

  CPPUNIT_TEST_SUITE_END();

  void testManualCreation1();
  void testManualCreation2();
  void testManualCreation3();
  void testPlatformException();

  //void testFail();

public:
  void setUp();
  void tearDown();
};

#endif // PLATFORMTEST_H

