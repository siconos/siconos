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
#include "LinearTIRXMLTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(LinearTIRXMLTest);




void LinearTIRXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("LinearTIR.xml");
    this->root = xmlDocGetRootElement(doc);
    this->LinearTIR = LinearTIRXML(root);
  }
  catch (SiconosException e)
  {
    cout << "Error in LinearTIRXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void LinearTIRXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void LinearTIRXMLTest::testGetC()
{
  try
  {
    SiconosMatrix m;
    m = LinearTIR.getC();

    CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : m.size(0)", m.size(0) == 2);
    CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : m.size(0)", m.size(1) == 3);
    CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : m", m(0, 2) == 3);
    cout << "LinearTIRXMLTest >>> testGetC ..................................... OK\n ";
  }
  catch (SiconosException e)
  {
    cout << "LinearTIRXMLTest >>> testGetC .................................... failed!\n" << e.report() << endl;
    exit(0);
  }
}

void LinearTIRXMLTest::testGetD()
{
  SiconosMatrix m;
  m = LinearTIR.getD();

  CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : m.size(0)", m.size(0) == 2);
  CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : m.size(0)", m.size(1) == 3);
  CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : m", m(0, 2) == 3);
  cout << "LinearTIRXMLTest >>> testGetD ..................................... OK\n ";
}

void LinearTIRXMLTest::testGeta()
{
  /*SiconosVector*/SimpleVector a;
  a = LinearTIR.getA();

  CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : a.size()", a.size() == 8);
  CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : a(0)", a(0) == 0.0);
  CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : a(1)", a(1) == 5.555);
  CPPUNIT_ASSERT_MESSAGE("testGetSiconosMatrixValue : a(2)", a(2) == -5.555);
  cout << "LinearTIRXMLTest >>> testGetA ..................................... OK\n ";
}

void LinearTIRXMLTest::testIfAttributeNotPresent()
{
  SiconosMatrix m;
  m = LinearTIR.getF();
}


void LinearTIRXMLTest::testIfTagIsNotPresent()
{
  xmlFreeDoc(doc);
  xmlCleanupParser();
  this->doc = xmlParseFile("LinearTIRCTagIsNotPresent.xml");
  this->root = xmlDocGetRootElement(doc);
  this->LinearTIR = LinearTIRXML(root);


}


