/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "LagrangianLinearTIRTest.hpp"
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
  e.reset(new SimpleVector(1));
  (*e)(0) = 0.1;
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
  cout << "--> Constructor xml test ended with success." << endl;
}

// data constructor (1)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR1()
{
  cout << "--> Test: constructor 1." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1c : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1d : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 1 test ended with success." << endl;
}

// data constructor (5)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR2()
{
  cout << "--> Test: constructor 2." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C, D, F, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2g : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 2 test ended with success." << endl;
}

// data constructor (5)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR3()
{
  cout << "--> Test: constructor 3." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3d : ", folr->e() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3g : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 3 test ended with success." << endl;
}

// data constructor (4)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR4()
{
  cout << "--> Test: constructor 4." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR4c : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR4d : ", folr->getSubType() == RELATION::LinearTIR, true);
  cout << "--> Constructor 4 test ended with success." << endl;
}

// // data constructor (5)
// void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR5()
// {
//   cout << "--> Test: constructor 5: obsolet" << endl;
// //   SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C,*D,*F,*e));
// //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5f : ", folr->getType()==RELATION::Lagrangian, true);
// //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR5g : ", folr->getSubType()==RELATION::LinearTIR, true);
//    cout << "--> Constructor 5 test ended with success." << endl;
// }

// // data constructor (6)
// void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR6()
// {
//   cout << "--> Test: constructor 6: obsolet" << endl;
// //   SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C,*e));
// //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR6f : ", folr->getType()==RELATION::Lagrangian, true);
// //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR6g : ", folr->getSubType()==RELATION::LinearTIR, true);
//   cout << "--> Constructor 5 test ended with success." << endl;
// }

// set C as a matrix and then plug it


// setCPtr
void LagrangianLinearTIRTest::testSetCPtr()
{
  cout << "--> Test: setCPtr." << endl;
  SP::SiconosMatrix tmp(new SimpleMatrix(*C));
  tmp->zero();
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*tmp));
  folr->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C() == C, true);
  cout << "--> setCPtr test ended with success." << endl;
}

// set D

// setDPtr
void LagrangianLinearTIRTest::testSetDPtr()
{
  cout << "--> Test: setDPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->D() == D, true);
  cout << "--> setDPtr test ended with success." << endl;
}

// set F

// setFPtr
void LagrangianLinearTIRTest::testSetFPtr()
{
  cout << "--> Test: setFPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->F() == F, true);
  cout << "--> setFPtr test ended with success." << endl;
}

// set E

// setEPtr
void LagrangianLinearTIRTest::testSetEPtr()
{
  cout << "--> Test: setEPtr." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(*C));
  folr->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->e() == e, true);
  cout << "--> setEPtr test ended with success." << endl;
}



void LagrangianLinearTIRTest::testGetJacPtr()
{
  cout << "--> Test: jac." << endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJachq: ", folr->jachq() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJachlambda: ", folr->jachlambda() == D, true);

  cout << "--> setBPtr test ended with success." << endl;
}

void LagrangianLinearTIRTest::End()
{
  cout << "===========================================" << endl;
  cout << " ===== End of LagrangianLinearTIR Tests ===== " << endl;
  cout << "=========================================== " << endl;
}
