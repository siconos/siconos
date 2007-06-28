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
#include "LagrangianLinearRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearRTest);


void LagrangianLinearRTest::setUp()
{
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("LagrangianLinearR_test.xml");
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

  // look for NSDS node
  xmlNode * nodetmp = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Content");
  // get relation
  xmlNode * node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LagrangianLinearRelation");
  tmpxml1 = new LagrangianLinearRXML(node1);

  H = new SimpleMatrix("matH.dat", true);
  b = new SimpleVector(1);
  D = new SimpleMatrix(1, 1);
  F = new SimpleMatrix(1, 1);
  (*b)(0) = 12;
  (*D)(0, 0) = 13;
  (*F)(0, 0) = 8;
}

void LagrangianLinearRTest::tearDown()
{
  delete tmpxml1;
  delete b;
  delete H;
  delete D;
  delete F;
}

// xml constructor (1)
void LagrangianLinearRTest::testBuildLagrangianLinearR0()
{
  cout << "========================================" << endl;
  cout << "=== LagrangianLinearR tests start ...=== " << endl;
  cout << "========================================" << endl;
  LagrangianLinearR * llr = new LagrangianLinearR(tmpxml1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR0a : ", llr->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR0a : ", llr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR0b : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR0c : ", llr->getB() == *b, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR0d : ", llr->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR0e : ", llr->getF() == *F, true);

  delete llr;
  cout << " xml Constructor LLR 1 ok" << endl;
}

// data constructor (1)
void LagrangianLinearRTest::testBuildLagrangianLinearR1()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR1a : ", llr->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR1b : ", llr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearLagrangianR1c : ", llr->getH() == *H, true);
  delete llr;
  cout << " xml Constructor LLR 2 ok" << endl;
}

// data constructor (2)
void LagrangianLinearRTest::testBuildLagrangianLinearR2()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H, *b);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR2a : ", llr->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR2b : ", llr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearLagrangianR2c : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR2d : ", llr->getB() == *b, true);
  delete llr;
  cout << " Constructor LLR 2 ok" << endl;
}

// data constructor (3)
void LagrangianLinearRTest::testBuildLagrangianLinearR3()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H, *b, *D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR3a : ", llr->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR3b : ", llr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearLagrangianR3c : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR3d : ", llr->getB() == *b, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR3e : ", llr->getD() == *D, true);
  delete llr;
  cout << " Constructor LLR 3 ok" << endl;
}

// data constructor (4)
void LagrangianLinearRTest::testBuildLagrangianLinearR4()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H, *b, *D, *F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR4a : ", llr->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR4b : ", llr->getSubType() == "LinearR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearLagrangianR4c : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR4d : ", llr->getB() == *b, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR4e : ", llr->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR4f : ", llr->getF() == *F, true);
  delete llr;
  cout << " Constructor LLR 4 ok" << endl;
}

// set H
void LagrangianLinearRTest::testSetH()
{
  SiconosMatrix * tmp = new SimpleMatrix(*H);
  tmp->zero();
  LagrangianLinearR * llr = new LagrangianLinearR(*tmp, *b);
  llr->setH(*H);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetH : ", llr->getH() == *H, true);
  delete llr;
  delete tmp;
  cout << " testSetH ok" << endl;
}

// setHPtr
void LagrangianLinearRTest::testSetHPtr()
{
  SiconosMatrix * tmp = new SimpleMatrix(*H);
  tmp->zero();
  LagrangianLinearR * llr = new LagrangianLinearR(*tmp, *b);
  llr->setHPtr(H);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetHPtr : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetHPtr : ", llr->getHPtr() == H, true);
  delete llr;
  delete tmp;
  cout << " test setHPtr ok" << endl;
}

// setB
void LagrangianLinearRTest::testSetB()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  llr->setB(*b);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetb: ", llr->getB() == *b, true);
  delete llr;
  cout << " test setB ok" << endl;
}

// setBPtr
void LagrangianLinearRTest::testSetBPtr()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  llr->setBPtr(b);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", llr->getB() == *b, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", llr->getBPtr() == b, true);
  delete llr;
  cout << " test setBPtr ok" << endl;
}

// setD
void LagrangianLinearRTest::testSetD()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  llr->setD(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetb: ", llr->getD() == *D, true);
  delete llr;
  cout << " test setD ok" << endl;
}

// setDPtr
void LagrangianLinearRTest::testSetDPtr()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  llr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", llr->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", llr->getDPtr() == D, true);
  delete llr;
  cout << " test setDPtr ok" << endl;
}

void LagrangianLinearRTest::End()
{
  cout << "============================================" << endl;
  cout << " ===== End of LagrangianLinearR tests ===== " << endl;
  cout << "============================================" << endl;
}
