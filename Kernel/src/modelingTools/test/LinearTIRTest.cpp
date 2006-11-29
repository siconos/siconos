/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
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
#include "LinearTIRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LinearTIRTest);


void LinearTIRTest::setUp()
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LinearTimeInvariantRelation");
  tmpxml1 = new LinearTIRXML(node1);

  // second file
  // parse xml file:
  xmlDocPtr doc2;
  xmlNodePtr cur2;
  doc2 = xmlParseFile("linearTIR_test2.xml");
  if (doc2 == NULL)
    XMLException::selfThrow("Document not parsed successfully");
  cur2 = xmlDocGetRootElement(doc2);
  if (cur2 == NULL)
  {
    XMLException::selfThrow("empty document");
    xmlFreeDoc(doc2);
  }
  // get rootNode
  if (xmlStrcmp(cur2->name, (const xmlChar *) "SiconosModel"))
  {
    XMLException::selfThrow("document of the wrong type, root node !=SiconosModel");
    xmlFreeDoc(doc2);
  }
  // look for NSDS, Interaction and relation node
  xmlNode* nodetmp2 = SiconosDOMTreeTools::findNodeChild(cur2, "NSDS");
  nodetmp2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "Interaction_Definition");
  nodetmp2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "Interaction");
  nodetmp2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "Interaction_Content");
  // get relation
  node2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "LinearTimeInvariantRelation");
  tmpxml2 = new LinearTIRXML(node2);
}

void LinearTIRTest::tearDown()
{
  delete nsds;
  delete tmpxml1;
  delete tmpxml2;
  delete e;
  delete F;
  delete B;
  delete C;
  delete D;
}

// xml constructor
void LinearTIRTest::testBuildLinearTIR0()
{
  cout << "================================" << endl;
  cout << "=== LinearTIR tests start ...=== " << endl;
  cout << "================================" << endl;
  LinearTIR * ltir = new LinearTIR(tmpxml1);
  cout << "--> Test: constructor xml." << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0a : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0b : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0c : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0d : ", ltir->getE() == *e, true);
  delete ltir;
  cout << "--> Constructor xml test ended with success." << endl;
}

// data constructor (1)
void LinearTIRTest::testBuildLinearTIR1()
{
  cout << "--> Test: constructor 1." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR : ", ltir->getB() == *B, true);
  delete ltir;
  cout << "--> Constructor 1 test ended with success." << endl;
}

// data constructor (2)
void LinearTIRTest::testBuildLinearTIR2()
{
  cout << "--> Test: constructor 2." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *D, *F, *e, *B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRC : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRD : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRF : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRE : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRB : ", ltir->getB() == *B, true);
  delete ltir;
  cout << "--> Constructor 2 test ended with success." << endl;
}

// copy constructor
void LinearTIRTest::testBuildLinearTIR3()
{
  cout << "--> Test: constructor 3." << endl;
  Relation * ref = new LinearTIR(tmpxml1);
  LinearTIR * ltir = new LinearTIR(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3C : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3D : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3F : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3E : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3B : ", ltir->getB() == *B, true);
  delete ltir;
  delete ref;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// xml constructor with input/output as plug-in
void LinearTIRTest::testBuildLinearTIR4()
{
  cout << "--> Test: constructor 4." << endl;
  LinearTIR * ltir = new LinearTIR(tmpxml2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4b : ", ltir->getType() == "LinearTIR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4c : ", ltir->getComputeOutputName() == "TestPlugin:y", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4d : ", ltir->getComputeInputName() == "TestPlugin:R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4e : ", ltir->getCPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4f : ", ltir->getBPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4g : ", ltir->getDPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4h : ", ltir->getEPtr() == NULL, true);
  delete ltir;
  cout << "--> Constructor 4 test ended with success." << endl;
}


// set C
void LinearTIRTest::testSetC()
{
  cout << "--> Test: setC." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*C);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*tmp, *B);
  ltir->setC(*C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", ltir->getC() == *C, true);
  delete ltir;
  delete tmp;
  cout << "--> setC test ended with success." << endl;
}

// setCPtr
void LinearTIRTest::testSetCPtr()
{
  cout << "--> Test: setCPtr." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*C);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*tmp, *B);
  ltir->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getCPtr() == C, true);
  delete ltir;
  delete tmp;
  cout << "--> setCPtr test ended with success." << endl;
}

// set D
void LinearTIRTest::testSetD()
{
  cout << "--> Test: setD." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setD(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetD: ", ltir->getD() == *D, true);
  delete ltir;
  cout << "--> setD test ended with success." << endl;
}

// setDPtr
void LinearTIRTest::testSetDPtr()
{
  cout << "--> Test: setDPtr." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", ltir->getDPtr() == D, true);
  delete ltir;
  cout << "--> setDPtr test ended with success." << endl;
}

// set F
void LinearTIRTest::testSetF()
{
  cout << "--> Test: setF." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setF(*F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetF: ", ltir->getF() == *F, true);
  delete ltir;
  cout << "--> setF test ended with success." << endl;
}

// setFPtr
void LinearTIRTest::testSetFPtr()
{
  cout << "--> Test: setFPtr." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", ltir->getFPtr() == F, true);
  delete ltir;
  cout << "--> setFPtr test ended with success." << endl;
}

// set E
void LinearTIRTest::testSetE()
{
  cout << "--> Test: setE." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setE(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetE: ", ltir->getE() == *e, true);
  delete ltir;
  cout << "--> setE test ended with success." << endl;
}

// setEPtr
void LinearTIRTest::testSetEPtr()
{
  cout << "--> Test: setEPtr." << endl;
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", ltir->getEPtr() == e, true);
  delete ltir;
  cout << "--> setEPtr test ended with success." << endl;
}

// set B
void LinearTIRTest::testSetB()
{
  cout << "--> Test: setB." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*B);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*C, *tmp);
  ltir->setB(*B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB: ", ltir->getB() == *B, true);
  delete ltir;
  delete tmp;
  cout << "--> setB test ended with success." << endl;
}

// setBPtr
void LinearTIRTest::testSetBPtr()
{
  cout << "--> Test: setBPtr." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*B);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*C, *tmp);
  ltir->setBPtr(B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", ltir->getBPtr() == B, true);
  delete ltir;
  delete tmp;
  cout << "--> setBPtr test ended with success." << endl;
}

void LinearTIRTest::End()
{
  cout << "====================================" << endl;
  cout << " ===== End of LinearTIR Tests ===== " << endl;
  cout << "====================================" << endl;
}
