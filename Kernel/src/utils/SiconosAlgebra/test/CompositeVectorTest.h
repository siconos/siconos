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
#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosVector.h"
#include "SimpleVector.h"
#include "CompositeVector.h"
#include "SiconosMatrix.h"

#include <math.h>

using namespace std;

class CompositeVectorTest : public CppUnit::TestFixture
{


private:

  // test suite
  CPPUNIT_TEST_SUITE(CompositeVectorTest);

  // test list
  CPPUNIT_TEST(testBuildCompositeVector);
  CPPUNIT_TEST(testBuildCompositeVector1);
  CPPUNIT_TEST(testOperatorAccess);
  CPPUNIT_TEST(testSetValue);
  CPPUNIT_TEST(testGetValue);
  CPPUNIT_TEST(testGetValues);
  CPPUNIT_TEST(testSetValues);
  CPPUNIT_TEST(testAdd);
  CPPUNIT_TEST(testSize);
  CPPUNIT_TEST(testReadWrite);
  CPPUNIT_TEST(testOperatorPlusEqual);
  CPPUNIT_TEST(testOperatorEqual);
  CPPUNIT_TEST(testOperatorComp);
  CPPUNIT_TEST(testOperatorMultDivEqual);
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testSubtraction);
  CPPUNIT_TEST(testExternalOperatorPlusMoins);
  CPPUNIT_TEST(testExternalOperatorMultDiv);
  CPPUNIT_TEST(testExternalOperatorMultMat);
  CPPUNIT_TEST(testExternalOperatorMatTransMult);

  CPPUNIT_TEST_SUITE_END();

  void testBuildCompositeVector();
  void testBuildCompositeVector1();
  void testOperatorAccess();
  void testSetValue();
  void testGetValue();
  void testSetValues();
  void testGetValues();
  void testAdd();
  void testSize();
  void testReadWrite();

  void testOperatorPlusEqual();
  void testOperatorEqual();
  void testOperatorComp();
  void testOperatorMultDivEqual();
  void testAddition();
  void testSubtraction();
  void testExternalOperatorPlusMoins();
  void testExternalOperatorMultDiv();
  void testExternalOperatorMultMat();
  void testExternalOperatorMatTransMult();
  // \todo exception test

  SimpleVector * q;
  SiconosVector * compVect, *r;
  SimpleVector * simpleVect;
  std::vector<double> vq;
  std::vector<double> vdotq;
  CompositeVector * CV, *tmp;


public:
  void setUp();
  void tearDown();

};

#endif



