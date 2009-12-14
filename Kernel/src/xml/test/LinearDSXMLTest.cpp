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
#include "LinearDSXMLTest.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(LinearsystemDSXMLTest);




void LinearsystemDSXMLTest::setUp()
{
  doc = xmlParseFile("LinearDS.xml");
  root = xmlDocGetRootElement(doc);
  child = NULL;
  matrixRef = SiconosMatrix("matrix.dat", true);
  vectorRef = /*SiconosVector*/SimpleVector("vector.dat", true);
}

void LinearsystemDSXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}

void LinearsystemDSXMLTest::testGetAB()
{
  //  try{
  LinearDSXML lsdsxml(root, false);
  SiconosMatrix A, E;
  A = lsdsxml.getA();
  E = lsdsxml.getE();

  CPPUNIT_ASSERT_MESSAGE("testGetAB : A == matrixRef", A == matrixRef);
  CPPUNIT_ASSERT_MESSAGE("testGetAB : E == matrixRef", E == matrixRef);
  //  }
  //  catch(SiconosException e)
  //  {
  //    cout << "LinearsystemDSXMLTest error : "<<e.report() <<endl;
  //    exit(0);
  //  }

  cout << " LinearsystemDSXMLTest >>> testGetAB ................................ OK\n ";
}

void LinearsystemDSXMLTest::testGetUF()
{
  LinearDSXML lsdsxml(root, false);
  /*SiconosVector*/
  SimpleVector U, b;
  U = lsdsxml.getUVector();
  b = lsdsxml.getBVector();
  CPPUNIT_ASSERT_MESSAGE("testGetUF : U == vectorRef", U == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetUF : b == vectorRef", b == vectorRef);
  cout << " LinearsystemDSXMLTest >>> testGetUF ................................ OK\n ";
}
