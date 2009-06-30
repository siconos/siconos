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
#include "FirstOrderLinearRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearRTest);


void FirstOrderLinearRTest::setUp()
{
  C.reset(new SimpleMatrix("matC.dat", true));
  D.reset(new SimpleMatrix("matD.dat", true));
  B.reset(new SimpleMatrix("matB.dat", true));
  F.reset(new SimpleMatrix("matF.dat", true));
  e.reset(new SimpleVector(2));
  Cp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:C"));
  Dp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:D"));
  Bp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:B"));
  Fp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:F"));
  ep.reset(new Plugged_Vector_FTime("TestPlugin:e"));
  (*e)(0) = 0.1;
  (*e)(1) = 0.1;
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("firstOrderLinearR_test.xml");
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderRelation");
  tmpxml1.reset(new LinearRXML(node1));

}

void FirstOrderLinearRTest::tearDown()
{}

// xml constructor
void FirstOrderLinearRTest::testBuildFirstOrderLinearR0()
{
  cout << "==========================================" << endl;
  cout << "==== FirstOrderLinearR tests start ...====" << endl;
  cout << "==========================================" << endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR(tmpxml1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getSubType() == RELATION::LinearR, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getCPtr()->getPluginName()=="TestPlugin:C", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getDPtr()->getPluginName()=="TestPlugin:D", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getFPtr()->getPluginName()=="TestPlugin:F", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getEPtr()->getPluginName()=="TestPlugin:e", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR0 : ", folr->getBPtr()->getPluginName()=="TestPlugin:B", true);
  cout << "--> Constructor xml test ended with success." << endl;
}

// data constructor (1)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR1()
{
  cout << "--> Test: constructor 1." << endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR("TestPlugin:C", "TestPlugin:B"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getSubType() == RELATION::LinearR, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getCPtr()->getPluginName()=="TestPlugin:C", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getBPtr()->getPluginName()=="TestPlugin:B", true);
  cout << "--> Constructor 1 test ended with success." << endl;
}

void FirstOrderLinearRTest::testBuildFirstOrderLinearR3()
{
  cout << "--> Test: constructor 3." << endl;

  SP::FirstOrderLinearR folr(new FirstOrderLinearR("TestPlugin:C", "TestPlugin:D", "TestPlugin:F", "TestPlugin:e", "TestPlugin:B"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getSubType() == RELATION::LinearR, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getCPtr()->getPluginName()=="TestPlugin:C", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getDPtr()->getPluginName()=="TestPlugin:D", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getFPtr()->getPluginName()=="TestPlugin:F", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getEPtr()->getPluginName()=="TestPlugin:e", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getBPtr()->getPluginName()=="TestPlugin:B", true);
  cout << "--> Constructor 3 test ended with success." << endl;
}

// data constructor (4)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR4()
{
  cout << "--> Test: constructor 4." << endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR(Cp, Bp));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4a : ", folr->getCPtr() == Cp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4b : ", folr->getBPtr() == Bp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4c : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4d : ", folr->getSubType() == RELATION::LinearR, true);
  cout << "--> Constructor 4 test ended with success." << endl;
}

// data constructor (5)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR5()
{
  cout << "--> Test: constructor 5." << endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR(Cp, Dp, Fp, ep, Bp));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5a : ", folr->getCPtr() == Cp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5b : ", folr->getDPtr() == Dp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5c : ", folr->getFPtr() == Fp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5d : ", folr->getEPtr() == ep, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5e : ", folr->getBPtr() == Bp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5f : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5g : ", folr->getSubType() == RELATION::LinearR, true);
  cout << "--> Constructor 5 test ended with success." << endl;
}

// set C as a matrix and then plug it
void FirstOrderLinearRTest::testSetC()
{
  //   cout << "--> Test: setC." << endl;
  //   SP::SiconosMatrix tmp(new SimpleMatrix(*C));
  //   tmp->zero();
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*tmp,*B));
  //   folr->setC(*C);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getC()==*C, true);
  //   //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getCPtr()->isPlugged()==false, true);

  //   folr->setComputeCFunction("TestPlugin.so","C");
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getCPtr()->isPlugged()==true, true);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getCPtr()->getPluginName()=="TestPlugin:C", true);
  //   cout << "--> setC test ended with success." << endl;
}

// setCPtr
void FirstOrderLinearRTest::testSetCPtr()
{
  //   cout << "--> Test: setCPtr." << endl;
  //   FirstOrderLinearR::SP_PluggedMatrix tmp(new FirstOrderLinearR::PluggedMatrix(*C));
  //   tmp->zero();
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*tmp,*B));
  //   folr->setCPtr(Cp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getC()==*Cp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getCPtr()==Cp, true);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getCPtr()->isPlugged(), true);
  // //   folr->setComputeCFunction("TestPlugin.so","C");
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getCPtr()->isPlugged()==true, true);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->getCPtr()->getPluginName()=="TestPlugin:C", true);

  cout << "--> setCPtr test ended with success." << endl;
}
// set C as a plugin
void FirstOrderLinearRTest::testSetCPtr2()
{
  //   cout << "--> Test: setCPtr2." << endl;
  //   FirstOrderLinearR::SP_PluggedMatrix tmp(new FirstOrderLinearR::PluggedMatrix("TestPlugin:D"));
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(tmp,Bp));
  //   folr->setCPtr(Cp);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC2 : ", folr->getCPtr()->isPlugged()==true, true);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC2 : ", folr->getCPtr()->getPluginName()=="TestPlugin:C", true);

  //   cout << "--> setCPtr2 test ended with success." << endl;
}

// set D
void FirstOrderLinearRTest::testSetD()
{
  //   cout << "--> Test: setD." << endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setD(*D);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetD: ", folr->getD()==*D, true);
  //   cout << "--> setD test ended with success." << endl;
}

// setDPtr
void FirstOrderLinearRTest::testSetDPtr()
{
  //   cout << "--> Test: setDPtr." << endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setDPtr(Dp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", folr->getD()==*Dp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->getDPtr()==Dp, true);
  //   cout << "--> setDPtr test ended with success." << endl;
}

// set F
void FirstOrderLinearRTest::testSetF()
{
  //   cout << "--> Test: setF." << endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setF(*F);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetF: ", folr->getF()==*F, true);
  //   cout << "--> setF test ended with success." << endl;
}

// setFPtr
void FirstOrderLinearRTest::testSetFPtr()
{
  //   cout << "--> Test: setFPtr." << endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setFPtr(Fp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", folr->getF()==*Fp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->getFPtr()==Fp, true);
  //   cout << "--> setFPtr test ended with success." << endl;
}

// set E
void FirstOrderLinearRTest::testSetE()
{
  cout << "--> Test: setE." << endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C, *B));
  folr->setE(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetE: ", folr->getE() == *e, true);
  cout << "--> setE test ended with success." << endl;
}

// setEPtr
void FirstOrderLinearRTest::testSetEPtr()
{
  //   cout << "--> Test: setEPtr." << endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setEPtr(ep);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", folr->getE()==*ep, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->getEPtr()==ep, true);
  //   cout << "--> setEPtr test ended with success." << endl;
}

// set B
void FirstOrderLinearRTest::testSetB()
{
  //   cout << "--> Test: setB." << endl;
  //   SP::SiconosMatrix tmp(new SimpleMatrix(*B));
  //   tmp->zero();
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*tmp));
  //   folr->setB(*B);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB: ", folr->getB()==*B, true);
  //   cout << "--> setB test ended with success." << endl;
}

// setBPtr
void FirstOrderLinearRTest::testSetBPtr()
{
  //   cout << "--> Test: setBPtr." << endl;
  //   FirstOrderLinearR::SP_PluggedMatrix tmp(new FirstOrderLinearR::PluggedMatrix(*B));
  //   tmp->zero();
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*tmp));
  //   folr->setBPtr(Bp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", folr->getB()==*Bp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", folr->getBPtr()==Bp, true);
  //   cout << "--> setBPtr test ended with success." << endl;
}

void FirstOrderLinearRTest::End()
{
  cout << "===========================================" << endl;
  cout << " ===== End of FirstOrderLinearR Tests ===== " << endl;
  cout << "=========================================== " << endl;
}
