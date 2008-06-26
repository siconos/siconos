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
#include "FirstOrderLinearTIRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearTIRTest);


void FirstOrderLinearTIRTest::setUp()
{
  C = new SimpleMatrix("matC.dat", true);
  D = new SimpleMatrix("matD.dat", true);
  B = new SimpleMatrix("matB.dat", true);
  F = new SimpleMatrix("matF.dat", true);
  e = new SimpleVector(2);
  (*e)(0) = 0.1;
  (*e)(1) = 0.1;
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("firstOrderLinearTIR_test.xml");
  if (doc == NULL)
    XMLException::selfThrow("Document not parsed successfully");
  cur = xmlDocGetRootElement(doc);
  if (cur == NULL)
  {
    XMLException::selfThrow("empty document");
    xmlFreeDoc(doc);
  }

  // get rootNode

  if (xmlStrcmp(cur->name, (const xmlChar *) "SiconosModel"))
  {
    XMLException::selfThrow("document of the wrong type, root node !=SiconosModel");
    xmlFreeDoc(doc);
  }

  // look for NSDS, Interaction and relation node
  xmlNode* nodetmp = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  NonSmoothDynamicalSystemXML * nsdsxml = new NonSmoothDynamicalSystemXML(nodetmp);
  nsds = new NonSmoothDynamicalSystem(nsdsxml);
  delete nsdsxml;
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Definition");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Content");
  // get relation
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderRelation");
  tmpxml1 = new FirstOrderLinearRXML(node1);
}

void FirstOrderLinearTIRTest::tearDown()
{
  delete nsds;
  delete tmpxml1;
  delete e;
  delete F;
  delete B;
  delete C;
  delete D;
}

// xml constructor
void FirstOrderLinearTIRTest::testBuildFirstOrderLinearTIR0()
{
  cout << "================================" << endl;
  cout << "=== FirstOrderLinearTIR tests start ...=== " << endl;
  cout << "================================" << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(tmpxml1);
  cout << "--> Test: constructor xml." << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0a : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0b : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0c : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0d : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0e : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0f : ", ltir->getType() == FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0g : ", ltir->getSubType() == LinearTIR, true);
  delete ltir;
  cout << "--> Constructor xml test ended with success." << endl;
}

// data constructor (1)
void FirstOrderLinearTIRTest::testBuildFirstOrderLinearTIR1()
{
  cout << "--> Test: constructor 1." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1a : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1b : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1c : ", ltir->getType() == FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1d : ", ltir->getSubType() == LinearTIR, true);
  delete ltir;
  cout << "--> Constructor 1 test ended with success." << endl;
}

// data constructor (2)
void FirstOrderLinearTIRTest::testBuildFirstOrderLinearTIR2()
{
  cout << "--> Test: constructor 2." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *D, *F, *e, *B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2a : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2b : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2c : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2d : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2e : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2f : ", ltir->getType() == FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2g : ", ltir->getSubType() == LinearTIR, true);
  delete ltir;
  cout << "--> Constructor 2 test ended with success." << endl;
}

void FirstOrderLinearTIRTest::End()
{
  cout << "====================================" << endl;
  cout << " ===== End of FirstOrderLinearTIR Tests ===== " << endl;
  cout << "====================================" << endl;
}
