/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
  doc = xmlParseFile("linearTIR_test.xml");
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderLinearTimeInvariantRelation");
  tmpxml1 = new FirstOrderLinearTIRXML(node1);
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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0f : ", ltir->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR0g : ", ltir->getSubType() == "LinearTIR", true);
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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1c : ", ltir->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1d : ", ltir->getSubType() == "LinearTIR", true);
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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2f : ", ltir->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2g : ", ltir->getSubType() == "LinearTIR", true);
  delete ltir;
  cout << "--> Constructor 2 test ended with success." << endl;
}

// set C
void FirstOrderLinearTIRTest::testSetC()
{
  cout << "--> Test: setC." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*C);
  tmp->zero();
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*tmp, *B);
  ltir->setC(*C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", ltir->getC() == *C, true);
  delete ltir;
  delete tmp;
  cout << "--> setC test ended with success." << endl;
}

// setCPtr
void FirstOrderLinearTIRTest::testSetCPtr()
{
  cout << "--> Test: setCPtr." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*C);
  tmp->zero();
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*tmp, *B);
  ltir->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getCPtr() == C, true);
  delete ltir;
  delete tmp;
  cout << "--> setCPtr test ended with success." << endl;
}

// set D
void FirstOrderLinearTIRTest::testSetD()
{
  cout << "--> Test: setD." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *B);
  ltir->setD(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetD: ", ltir->getD() == *D, true);
  delete ltir;
  cout << "--> setD test ended with success." << endl;
}

// setDPtr
void FirstOrderLinearTIRTest::testSetDPtr()
{
  cout << "--> Test: setDPtr." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *B);
  ltir->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", ltir->getDPtr() == D, true);
  delete ltir;
  cout << "--> setDPtr test ended with success." << endl;
}

// set F
void FirstOrderLinearTIRTest::testSetF()
{
  cout << "--> Test: setF." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *B);
  ltir->setF(*F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetF: ", ltir->getF() == *F, true);
  delete ltir;
  cout << "--> setF test ended with success." << endl;
}

// setFPtr
void FirstOrderLinearTIRTest::testSetFPtr()
{
  cout << "--> Test: setFPtr." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *B);
  ltir->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", ltir->getFPtr() == F, true);
  delete ltir;
  cout << "--> setFPtr test ended with success." << endl;
}

// set E
void FirstOrderLinearTIRTest::testSetE()
{
  cout << "--> Test: setE." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *B);
  ltir->setE(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetE: ", ltir->getE() == *e, true);
  delete ltir;
  cout << "--> setE test ended with success." << endl;
}

// setEPtr
void FirstOrderLinearTIRTest::testSetEPtr()
{
  cout << "--> Test: setEPtr." << endl;
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *B);
  ltir->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", ltir->getEPtr() == e, true);
  delete ltir;
  cout << "--> setEPtr test ended with success." << endl;
}

// set B
void FirstOrderLinearTIRTest::testSetB()
{
  cout << "--> Test: setB." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*B);
  tmp->zero();
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *tmp);
  ltir->setB(*B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB: ", ltir->getB() == *B, true);
  delete ltir;
  delete tmp;
  cout << "--> setB test ended with success." << endl;
}

// setBPtr
void FirstOrderLinearTIRTest::testSetBPtr()
{
  cout << "--> Test: setBPtr." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*B);
  tmp->zero();
  FirstOrderLinearTIR * ltir = new FirstOrderLinearTIR(*C, *tmp);
  ltir->setBPtr(B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", ltir->getBPtr() == B, true);
  delete ltir;
  delete tmp;
  cout << "--> setBPtr test ended with success." << endl;
}

void FirstOrderLinearTIRTest::End()
{
  cout << "====================================" << endl;
  cout << " ===== End of FirstOrderLinearTIR Tests ===== " << endl;
  cout << "====================================" << endl;
}
