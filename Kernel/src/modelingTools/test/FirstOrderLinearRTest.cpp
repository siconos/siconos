/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#include "FirstOrderLinearRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearRTest);


void FirstOrderLinearRTest::setUp()
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
  doc = xmlParseFile("firstOrderLinearR_test.xml");
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderLinearRelation");
  tmpxml1 = new FirstOrderLinearRXML(node1);
}

void FirstOrderLinearRTest::tearDown()
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
void FirstOrderLinearRTest::testBuildFirstOrderLinearR0()
{
  cout << "==========================================" << endl;
  cout << "==== FirstOrderLinearR tests start ...====" << endl;
  cout << "==========================================" << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getFunctionName("C") == "TestPlugin:C", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getFunctionName("D") == "TestPlugin:D", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getFunctionName("F") == "TestPlugin:F", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getFunctionName("e") == "TestPlugin:e", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getFunctionName("B") == "TestPlugin:B", true);
  delete folr;
  cout << "--> Constructor xml test ended with success." << endl;
}

// data constructor (1)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR1()
{
  cout << "--> Test: constructor 1." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR("TestPlugin:C", "TestPlugin:B");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getFunctionName("C") == "TestPlugin:C", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getFunctionName("B") == "TestPlugin:B", true);
  delete folr;
  cout << "--> Constructor 1 test ended with success." << endl;
}

// data constructor (2)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR2()
{
  cout << "--> Test: constructor 2." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR("TestPlugin:e");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR2 : ", folr->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR2 : ", folr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR2 : ", folr->getFunctionName("e") == "TestPlugin:e", true);
  delete folr;
  cout << "--> Constructor 2 test ended with success." << endl;
}

void FirstOrderLinearRTest::testBuildFirstOrderLinearR3()
{
  cout << "--> Test: constructor 3." << endl;

  FirstOrderLinearR * folr = new FirstOrderLinearR("TestPlugin:C", "TestPlugin:D", "TestPlugin:F", "TestPlugin:e", "TestPlugin:B");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getFunctionName("C") == "TestPlugin:C", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getFunctionName("D") == "TestPlugin:D", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getFunctionName("F") == "TestPlugin:F", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getFunctionName("e") == "TestPlugin:e", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getFunctionName("B") == "TestPlugin:B", true);
  delete folr;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// data constructor (4)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR4()
{
  cout << "--> Test: constructor 4." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4a : ", folr->getCPtr() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4b : ", folr->getBPtr() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4c : ", folr->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4d : ", folr->getSubType() == "LinearR", true);
  delete folr;
  cout << "--> Constructor 4 test ended with success." << endl;
}

// data constructor (5)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR5()
{
  cout << "--> Test: constructor 5." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, D, F, e, B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5a : ", folr->getCPtr() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5b : ", folr->getDPtr() == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5c : ", folr->getFPtr() == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5d : ", folr->getEPtr() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5e : ", folr->getBPtr() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5f : ", folr->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5g : ", folr->getSubType() == "LinearR", true);
  delete folr;
  cout << "--> Constructor 5 test ended with success." << endl;
}

// set C
void FirstOrderLinearRTest::testSetC()
{
  cout << "--> Test: setC." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*C);
  tmp->zero();
  FirstOrderLinearR * folr = new FirstOrderLinearR(tmp, B);
  folr->setC(*C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getC() == *C, true);
  delete folr;
  delete tmp;
  cout << "--> setC test ended with success." << endl;
}

// setCPtr
void FirstOrderLinearRTest::testSetCPtr()
{
  cout << "--> Test: setCPtr." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*C);
  tmp->zero();
  FirstOrderLinearR * folr = new FirstOrderLinearR(tmp, B);
  folr->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getCPtr() == C, true);
  delete folr;
  delete tmp;
  cout << "--> setCPtr test ended with success." << endl;
}

// set D
void FirstOrderLinearRTest::testSetD()
{
  cout << "--> Test: setD." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, B);
  folr->setD(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetD: ", folr->getD() == *D, true);
  delete folr;
  cout << "--> setD test ended with success." << endl;
}

// setDPtr
void FirstOrderLinearRTest::testSetDPtr()
{
  cout << "--> Test: setDPtr." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, B);
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", folr->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->getDPtr() == D, true);
  delete folr;
  cout << "--> setDPtr test ended with success." << endl;
}

// set F
void FirstOrderLinearRTest::testSetF()
{
  cout << "--> Test: setF." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, B);
  folr->setF(*F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetF: ", folr->getF() == *F, true);
  delete folr;
  cout << "--> setF test ended with success." << endl;
}

// setFPtr
void FirstOrderLinearRTest::testSetFPtr()
{
  cout << "--> Test: setFPtr." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, B);
  folr->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", folr->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->getFPtr() == F, true);
  delete folr;
  cout << "--> setFPtr test ended with success." << endl;
}

// set E
void FirstOrderLinearRTest::testSetE()
{
  cout << "--> Test: setE." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, B);
  folr->setE(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetE: ", folr->getE() == *e, true);
  delete folr;
  cout << "--> setE test ended with success." << endl;
}

// setEPtr
void FirstOrderLinearRTest::testSetEPtr()
{
  cout << "--> Test: setEPtr." << endl;
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, B);
  folr->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", folr->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->getEPtr() == e, true);
  delete folr;
  cout << "--> setEPtr test ended with success." << endl;
}

// set B
void FirstOrderLinearRTest::testSetB()
{
  cout << "--> Test: setB." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*B);
  tmp->zero();
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, tmp);
  folr->setB(*B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB: ", folr->getB() == *B, true);
  delete folr;
  delete tmp;
  cout << "--> setB test ended with success." << endl;
}

// setBPtr
void FirstOrderLinearRTest::testSetBPtr()
{
  cout << "--> Test: setBPtr." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*B);
  tmp->zero();
  FirstOrderLinearR * folr = new FirstOrderLinearR(C, tmp);
  folr->setBPtr(B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", folr->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", folr->getBPtr() == B, true);
  delete folr;
  delete tmp;
  cout << "--> setBPtr test ended with success." << endl;
}

void FirstOrderLinearRTest::End()
{
  cout << "===========================================" << endl;
  cout << " ===== End of FirstOrderLinearR Tests ===== " << endl;
  cout << "=========================================== " << endl;
}
