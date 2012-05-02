/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#ifndef __OSNSMatrixTest__
#define __OSNSMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "OSNSMatrix.hpp"
#include "Model.hpp"
#include "TimeStepping.hpp"
#include "OneStepNSProblem.hpp"

class OSNSMatrixTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(OSNSMatrixTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(OSNSMatrixTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildOSNSMatrix0);
  CPPUNIT_TEST(testBuildOSNSMatrix1);
  CPPUNIT_TEST(testBuildOSNSMatrix2);
  CPPUNIT_TEST(testFill);
  CPPUNIT_TEST(testConvert);
  CPPUNIT_TEST(testFill2);
  CPPUNIT_TEST(testConvert2);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildOSNSMatrix0();
  void testBuildOSNSMatrix1();
  void testBuildOSNSMatrix2();
  void testFill();
  void testConvert();
  void testFill2();
  void testConvert2();
  void End();
  // Members

  unsigned int n;
  double tol;
  SP::InteractionsSet indexSet;
  MapOfMapOfInteractionMatrices blocks;
  Model * temp ;
public:
  void setUp();
  void tearDown();

};

#endif




