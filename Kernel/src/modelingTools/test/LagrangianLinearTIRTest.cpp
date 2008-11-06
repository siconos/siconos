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
#include "LagrangianLinearTIRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearTIRTest);


void LagrangianLinearTIRTest::setUp()
{
  C.reset(new SimpleMatrix("matC.dat", true));
  D.reset(new SimpleMatrix("matD.dat", true));
  F.reset(new SimpleMatrix("matF.dat", true));
  e.reset(new SimpleVector(2));
  (*e)(0) = 0.1;
  (*e)(1) = 0.1;
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("LagrangianLinearTIR_test.xml");
  if (!doc)
    XMLException::selfThrow("Document not parsed successfully");
  cur = xmlDocGetRootElement(doc);
  if (!cur)
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
  SP::NonSmoothDynamicalSystemXML nsdsxml(new NonSmoothDynamicalSystemXML(nodetmp));
  nsds.reset(new NonSmoothDynamicalSystem(nsdsxml));
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Definition");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Content");

  // get relation
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LagrangianRelation");
  tmpxml1.reset(new LinearRXML(node1));

}

void LagrangianLinearTIRTest::tearDown()
{}

// xml constructor
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR0()
{
  cout << "==========================================" << endl;
  cout << "==== LagrangianLinearTIR tests start ...====" << endl;
  cout << "==========================================" << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(tmpxml1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR0 : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR0 : ", folr->getSubType() == RELATION::LinearTIR, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR0 : ", folr->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR0 : ", folr->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR0 : ", folr->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR0 : ", folr->getE() == *e, true);
  cout << "--> Constructor xml test ended with success." << endl;
}

// data constructor (1)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR1()
{
  cout << "--> Test: constructor 1." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1a : ", folr->getCPtr() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1c : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1d : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 1 test ended with success." << endl;
}

// data constructor (5)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR2()
{
  cout << "--> Test: constructor 2." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C, D, F, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2a : ", folr->getCPtr() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2b : ", folr->getDPtr() == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2c : ", folr->getFPtr() == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2d : ", folr->getEPtr() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2g : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 2 test ended with success." << endl;
}

// data constructor (5)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR3()
{
  cout << "--> Test: constructor 3." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3a : ", folr->getCPtr() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3d : ", folr->getEPtr() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3g : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 3 test ended with success." << endl;
}

// data constructor (4)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR4()
{
  cout << "--> Test: constructor 4." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR4a : ", folr->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR4c : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR4d : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 4 test ended with success." << endl;
}

// data constructor (5)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR5()
{
  cout << "--> Test: constructor 5." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C, *D, *F, *e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5a : ", folr->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5b : ", folr->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5c : ", folr->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5d : ", folr->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5g : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 5 test ended with success." << endl;
}

// data constructor (6)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR6()
{
  cout << "--> Test: constructor 6." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C, *e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR6a : ", folr->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR6d : ", folr->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR6f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR6g : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 5 test ended with success." << endl;
}

// set C as a matrix and then plug it
void LagrangianLinearTIRTest::testSetC()
{
  cout << "--> Test: setC." << endl;
  SP::SiconosMatrix tmp(new SimpleMatrix(*C));
  tmp->zero();
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*tmp));
  folr->setC(*C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getC() == *C, true);
  cout << "--> setC test ended with success." << endl;
}

// setCPtr
void LagrangianLinearTIRTest::testSetCPtr()
{
  cout << "--> Test: setCPtr." << endl;
  SP::SiconosMatrix tmp(new SimpleMatrix(*C));
  tmp->zero();
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*tmp));
  folr->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getCPtr() == C, true);
  cout << "--> setCPtr test ended with success." << endl;
}

// set D
void LagrangianLinearTIRTest::testSetD()
{
  cout << "--> Test: setD." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setD(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetD: ", folr->getD() == *D, true);
  cout << "--> setD test ended with success." << endl;
}

// setDPtr
void LagrangianLinearTIRTest::testSetDPtr()
{
  cout << "--> Test: setDPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->getDPtr() == D, true);
  cout << "--> setDPtr test ended with success." << endl;
}

// set F
void LagrangianLinearTIRTest::testSetF()
{
  cout << "--> Test: setF." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setF(*F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetF: ", folr->getF() == *F, true);
  cout << "--> setF test ended with success." << endl;
}

// setFPtr
void LagrangianLinearTIRTest::testSetFPtr()
{
  cout << "--> Test: setFPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->getFPtr() == F, true);
  cout << "--> setFPtr test ended with success." << endl;
}

// set E
void LagrangianLinearTIRTest::testSetE()
{
  cout << "--> Test: setE." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setE(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetE: ", folr->getE() == *e, true);
  cout << "--> setE test ended with success." << endl;
}

// setEPtr
void LagrangianLinearTIRTest::testSetEPtr()
{
  cout << "--> Test: setEPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->getEPtr() == e, true);
  cout << "--> setEPtr test ended with success." << endl;
}


void LagrangianLinearTIRTest::testGetJac()
{
  cout << "--> Test: getJacPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJacH: ", folr->getJacH(0) == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJacH: ", folr->getJacH(1) == *D, true);

  cout << "--> setBPtr test ended with success." << endl;
}

void LagrangianLinearTIRTest::testGetJacPtr()
{
  cout << "--> Test: getJacPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJacH: ", folr->getJacHPtr(0) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJacH: ", folr->getJacHPtr(1) == D, true);

  cout << "--> setBPtr test ended with success." << endl;
}

void LagrangianLinearTIRTest::End()
{
  cout << "===========================================" << endl;
  cout << " ===== End of LagrangianLinearTIR Tests ===== " << endl;
  cout << "=========================================== " << endl;
}
