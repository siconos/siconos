/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
  C = new SiconosMatrix("matC.dat", true);
  D = new SiconosMatrix("matD.dat", true);
  B = new SiconosMatrix("matB.dat", true);
  F = new SiconosMatrix("matF.dat", true);
  a = new SimpleVector(3);
  (*a)(0) = 4;
  (*a)(1) = 5;
  (*a)(2) = 6;
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
  delete a;
  delete F;
  delete B;
  delete C;
  delete D;
}

// default constructor

// Default  -> private -> no test
/*void LinearTIRTest::testBuildLinearTIR()
{
  LinearTIR * ltir = new LinearTIR();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", ltir->getType()==LINEARTIRELATION, true);
  delete ltir;
  cout << " Constructor LTIR 0 ok" << endl;
}
*/
// xml constructor
void LinearTIRTest::testBuildLinearTIR0()
{
  cout << "========================================" << endl;
  cout << "=== LinearTIR tests start ...=== " << endl;
  cout << "========================================" << endl;
  LinearTIR * ltir = new LinearTIR(tmpxml1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0a : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0b : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0c : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0d : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR0e : ", ltir->getA() == *a, true);
  delete ltir;
  cout << " xml Constructor LTIR ok" << endl;
}

// data constructor (1)
void LinearTIRTest::testBuildLinearTIR1()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR : ", ltir->getB() == *B, true);
  delete ltir;
  cout << " Constructor LTIR 1 ok" << endl;

}

// data constructor (2)
void LinearTIRTest::testBuildLinearTIR2()
{
  LinearTIR * ltir = new LinearTIR(*C, *D, *F, *e, *B, *a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRC : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRD : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRF : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRE : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRB : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRA : ", ltir->getA() == *a, true);
  delete ltir;
  cout << " Constructor LTIR 2 ok" << endl;
}

// copy constructor
void LinearTIRTest::testBuildLinearTIR3()
{
  Relation * ref = new LinearTIR(tmpxml1);
  LinearTIR * ltir = new LinearTIR(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3C : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3D : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3F : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3E : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3B : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR3A : ", ltir->getA() == *a, true);
  delete ltir;
  delete ref;
  cout << " copy Constructor LTIR ok" << endl;
}

// xml constructor with input/output as plug-in
void LinearTIRTest::testBuildLinearTIR4()
{
  LinearTIR * ltir = new LinearTIR(tmpxml2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4b : ", ltir->getType() == "LinearTIR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4c : ", ltir->getComputeOutputName() == "TestPlugin:y", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4d : ", ltir->getComputeInputName() == "TestPlugin:R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4e : ", ltir->getCPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4f : ", ltir->getBPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4g : ", ltir->getDPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4h : ", ltir->getEPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR4i : ", ltir->getAPtr() == NULL, true);
  delete ltir;
  cout << " xml Constructor (2) LTIR ok" << endl;
}


// set C
void LinearTIRTest::testSetC()
{
  SiconosMatrix * tmp = new SiconosMatrix(*C);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*tmp, *B);
  ltir->setC(*C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", ltir->getC() == *C, true);
  delete ltir;
  delete tmp;
  cout << " testSetC ok" << endl;
}

// setCPtr
void LinearTIRTest::testSetCPtr()
{
  SiconosMatrix * tmp = new SiconosMatrix(*C);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*tmp, *B);
  ltir->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getCPtr() == C, true);
  delete ltir;
  delete tmp;
  cout << " test setCPtr ok" << endl;
}

// set D
void LinearTIRTest::testSetD()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setD(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetD: ", ltir->getD() == *D, true);
  delete ltir;
  cout << " test setD ok" << endl;
}

// setDPtr
void LinearTIRTest::testSetDPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", ltir->getDPtr() == D, true);
  delete ltir;
  cout << " test setDPtr ok" << endl;
}

// set F
void LinearTIRTest::testSetF()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setF(*F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetF: ", ltir->getF() == *F, true);
  delete ltir;
  cout << " test setF ok" << endl;
}

// setFPtr
void LinearTIRTest::testSetFPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", ltir->getFPtr() == F, true);
  delete ltir;
  cout << " test setFPtr ok" << endl;
}

// set E
void LinearTIRTest::testSetE()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setE(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetE: ", ltir->getE() == *e, true);
  delete ltir;
  cout << " test setE ok" << endl;
}

// setEPtr
void LinearTIRTest::testSetEPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", ltir->getEPtr() == e, true);
  delete ltir;
  cout << " test setEPtr ok" << endl;
}

// set B
void LinearTIRTest::testSetB()
{
  SiconosMatrix * tmp = new SiconosMatrix(*B);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*C, *tmp);
  ltir->setB(*B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB: ", ltir->getB() == *B, true);
  delete ltir;
  delete tmp;
  cout << " test setB ok" << endl;
}

// setBPtr
void LinearTIRTest::testSetBPtr()
{
  SiconosMatrix * tmp = new SiconosMatrix(*B);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*C, *tmp);
  ltir->setBPtr(B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", ltir->getBPtr() == B, true);
  delete ltir;
  delete tmp;
  cout << " test setBPtr ok" << endl;
}

// set A
void LinearTIRTest::testSetA()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setA(*a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetA: ", ltir->getA() == *a, true);
  delete ltir;
  cout << " test setA ok" << endl;
}

// setAPtr
void LinearTIRTest::testSetAPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setAPtr(a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr : ", ltir->getA() == *a, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr: ", ltir->getAPtr() == a, true);
  delete ltir;
  cout << " test setAPtr ok" << endl;
}

void LinearTIRTest::End()
{
  cout << "====================================" << endl;
  cout << " ===== End of LinearTIR Tests ===== " << endl;
  cout << "====================================" << endl;
}
