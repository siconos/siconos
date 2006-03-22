/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#include "DynamicalSystemTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(DynamicalSystemTest);


void DynamicalSystemTest::setUp()
{
  x0 = new SimpleVector(3);
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  u0 = new SimpleVector(2);
  (*u0)(0) = 4;
  (*u0)(1) = 5;

  T0 = new SimpleMatrix("matT0.dat", true);
  J0 = new SimpleMatrix("matJ0.dat", true);

  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("ds_test.xml");
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
  xmlNode* nodetmp = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "DS_Definition");
  // get first ds
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "DynamicalSystem");
  tmpxml1 = new DynamicalSystemXML(node1, false);
  // get second ds
  node2 = SiconosDOMTreeTools::findFollowNode(node1, "DynamicalSystem");
  tmpxml2 = new DynamicalSystemXML(node2, false);
}

void DynamicalSystemTest::tearDown()
{
  delete tmpxml1;
  delete tmpxml2;
  delete x0;
  delete u0;
  delete T0;
  delete J0;
}

// xml constructor (1), without plugin
void DynamicalSystemTest::testBuildDynamicalSystem1()
{
  cout << "======================================" << endl;
  cout << "=== DynamicalSystem tests start ...=== " << endl;
  cout << "======================================" << endl;
  cout << "--> Test: constructor xml." << endl;
  DynamicalSystem * ds = new DynamicalSystem(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1A : ", ds->getType() == NLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1C : ", ds->getId() == "testDS1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1D : ", ds->getStepsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1E : ", ds->getX0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1F : ", ds->getU() == *u0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1G : ", ds->getT() == *T0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1H : ", ds->getIsPlugin()[0] == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem1I : ", ds->getIsPlugin()[1] == false, true);

  delete ds;
  cout << "--> Constructor xml test ended with success." << endl;
}


// xml constructor (2), with plugins
void DynamicalSystemTest::testBuildDynamicalSystem2()
{
  cout << "--> Test: constructor xml 2." << endl;
  DynamicalSystem * ds = new DynamicalSystem(tmpxml2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2A : ", ds->getType() == NLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2B : ", ds->getNumber() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2C : ", ds->getId() == "testDS2", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2E : ", ds->getX0() == 2 * *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2F : ", ds->getIsPlugin()[0] == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2G : ", ds->getIsPlugin()[1] == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2H : ", ds->getUSize() == 2, true);

  double time = 1.5;
  ds->computeU(time);
  ds->computeT();
  SimpleVector * x01 = new SimpleVector(2);
  (*x01)(0) = (*x0)(0);
  (*x01)(1) = (*x0)(1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2I : ", ds->getU() == 2 * time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem2J : ", ds->getT() == *T0, true);
  delete x01;
  delete ds;
  cout << "--> Constructor xml 2 test ended with success." << endl;
}


// constructor from data
void DynamicalSystemTest::testBuildDynamicalSystem3()
{
  cout << "--> Test: constructor 3." << endl;
  DynamicalSystem * ds = new DynamicalSystem(13, 3, *x0, "TestPlugin:vectorField", "TestPlugin:computeJacobianX");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3A : ", ds->getType() == NLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3C : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3D : ", ds->getX0() == *x0, true);

  double time = 1.5;
  ds->computeVectorField(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3E : ", ds->getXDot() == time* *x0, true);
  delete ds;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// copy constructor
void DynamicalSystemTest::testBuildDynamicalSystem4()
{
  cout << "--> Test: constructor 4." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml1);
  DynamicalSystem * ds2 = new DynamicalSystem(*ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4A : ", ds2->getType() == NLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4B : ", ds2->getNumber() == -2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4C : ", ds2->getId() == "copy", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4D : ", ds2->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4E : ", ds2->getX0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4F : ", ds2->getU() == *u0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4G : ", ds2->getT() == *T0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4H : ", ds2->getIsPlugin()[0] == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem4I : ", ds2->getIsPlugin()[1] == false, true);

  delete ds1;
  delete ds2;
  cout << "--> Constructor 4 test ended with success." << endl;
}

// setX0
void DynamicalSystemTest::testSetX0()
{
  cout << "--> Test: setX0." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX0(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0 : ", ds1->getX0() == *x0, true);
  delete ds1;
  cout << "--> setX0 test ended with success." << endl;
}

// setX0Ptr
void DynamicalSystemTest::testSetX0Ptr()
{
  cout << "--> Test: setX0Ptr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX0Ptr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0Ptr : ", ds1->getX0() == *x0, true);
  delete ds1;
  cout << "--> setX0Ptr test ended with success." << endl;
}

// setX0 with exception
void DynamicalSystemTest::testSetX02()
{
  cout << "--> Test: setX02." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX0(*u0);
  cout << "--> setX02 test ended with success." << endl;
}

// setX
void DynamicalSystemTest::testSetX()
{
  cout << "--> Test: setX." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX : ", ds1->getX() == *x0, true);
  delete ds1;
  cout << "--> setX test ended with success." << endl;
}

// setXPtr
void DynamicalSystemTest::testSetXPtr()
{
  cout << "--> Test: setXPtr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXPtr : ", ds1->getX() == *x0, true);
  delete ds1;
  cout << "--> setXPtr test ended with success." << endl;
}

// setX with exception
void DynamicalSystemTest::testSetX2()
{
  cout << "--> Test: setX2." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX(*u0);
  delete ds1;
  cout << "--> setX2 test ended with success." << endl;
}

// setXDot
void DynamicalSystemTest::testSetXDot()
{
  cout << "--> Test: setXDot." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXDot(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXDot : ", ds1->getXDot() == *x0, true);
  delete ds1;
  cout << "--> setXDot test ended with success." << endl;
}

// setXDotPtr
void DynamicalSystemTest::testSetXDotPtr()
{
  cout << "--> Test: setXDotptr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXDotPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXDotPtr : ", ds1->getXDot() == *x0, true);
  delete ds1;
  cout << "--> setXDotPtr test ended with success." << endl;
}

// setXDot with exception
void DynamicalSystemTest::testSetXDot2()
{
  cout << "--> Test: setXDot2." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXDot(*u0);
  delete ds1;
  cout << "--> setXDot2 test ended with success." << endl;
}

// setXFree
void DynamicalSystemTest::testSetXFree()
{
  cout << "--> Test: setXFree." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXFree(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXFree : ", ds1->getXFree() == *x0, true);
  delete ds1;
  cout << "--> setXFree test ended with success." << endl;
}

// setXFreePtr
void DynamicalSystemTest::testSetXFreePtr()
{
  cout << "--> Test: setXFreePtr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXFreePtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXFreePtr : ", ds1->getXFree() == *x0, true);
  delete ds1;
  cout << "--> setXFreePtr test ended with success." << endl;
}

// setXFree with exception
void DynamicalSystemTest::testSetXFree2()
{
  cout << "--> Test: setXFree2." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXFree(*u0);
  delete ds1;
  cout << "--> setXFree2 test ended with success." << endl;
}

// setR
void DynamicalSystemTest::testSetR()
{
  cout << "--> Test: setR." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setR(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetR : ", ds1->getR() == *x0, true);
  delete ds1;
  cout << "--> setR test ended with success." << endl;
}

// setRPtr
void DynamicalSystemTest::testSetRPtr()
{
  cout << "--> Test: setRPtr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setRPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetRPtr : ", ds1->getR() == *x0, true);
  delete ds1;
  cout << "--> setRPtr test ended with success." << endl;
}

// setR with exception
void DynamicalSystemTest::testSetR2()
{
  cout << "--> Test: setR2." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setR(*u0);
  delete ds1;
  cout << "--> setR2 test ended with success." << endl;
}

// set JacobianX
void DynamicalSystemTest::testSetJacobianX()
{
  cout << "--> Test: setJacobianX." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setJacobianX(*J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianX : ", ds1->getJacobianX() == *J0, true);
  delete ds1;
  cout << "--> setJacobianX test ended with success." << endl;
}

// setJacobianXPtr
void DynamicalSystemTest::testSetJacobianXPtr()
{
  cout << "--> Test: setJacobianXPtr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setJacobianXPtr(J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianXPtr : ", ds1->getJacobianX() == *J0, true);
  delete ds1;
  cout << "--> setJacobianXPtr test ended with success." << endl;
}

// setJacobianX with exception
void DynamicalSystemTest::testSetJacobianX2()
{
  cout << "--> Test: setJacobianX2." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setJacobianX(*T0);
  delete ds1;
  cout << "--> setJacobianX2 test ended with success." << endl;
}

// setU
void DynamicalSystemTest::testSetU()
{
  cout << "--> Test: setU." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setU(*u0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetU : ", ds1->getU() == *u0, true);
  delete ds1;
  cout << "--> setU test ended with success." << endl;
}

// setUPtr
void DynamicalSystemTest::testSetUPtr()
{
  cout << "--> Test: setUPtr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setUPtr(u0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetUPtr : ", ds1->getU() == *u0, true);
  delete ds1;
  cout << "--> setUPtr test ended with success." << endl;
}

// setU with exception
void DynamicalSystemTest::testSetU2()
{
  cout << "--> Test: setU2." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setU(*x0);
  delete ds1;
  cout << "--> setU2 test ended with success." << endl;
}

// set T
void DynamicalSystemTest::testSetT()
{
  cout << "--> Test: setT." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setT(*T0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetT : ", ds1->getT() == *T0, true);
  delete ds1;
  cout << "--> setT test ended with success." << endl;
}

// setTPtr
void DynamicalSystemTest::testSetTPtr()
{
  cout << "--> Test: setTPtr." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setTPtr(T0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetTPtr : ", ds1->getT() == *T0, true);
  delete ds1;
  cout << "--> setTPtr test ended with success." << endl;
}

// setT with exception
void DynamicalSystemTest::testSetT2()
{
  cout << "--> Test: setT2." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setT(*J0);
  delete ds1;
  cout << "--> setT2 test ended with success." << endl;
}

// init
void DynamicalSystemTest::testInitMemory()
{
  cout << "--> Test: initMemory." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->initMemory(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem1 : ", ds1->getXMemoryPtr()->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem2 : ", ds1->getXDotMemoryPtr()->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem3 : ", ds1->getRMemoryPtr()->getMemorySize() == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem4 : ", ds1->getXMemoryPtr()->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem5 : ", ds1->getXDotMemoryPtr()->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem6 : ", ds1->getRMemoryPtr()->getNbVectorsInMemory() == 0, true);
  delete ds1;
  cout << "--> initMemory test ended with success." << endl;
}


// swap
void DynamicalSystemTest::testSwap()
{
  cout << "--> Test: swap." << endl;
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX(*x0);
  ds1->setXDot(*x0);
  ds1->setR(*x0);
  ds1->initMemory(1);
  ds1->swapInMemory();
  cout << " swap ok " << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap1 : ", *((ds1->getXMemoryPtr()->getVectorMemory())[0]) == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap2 : ", *((ds1->getXDotMemoryPtr()->getVectorMemory())[0]) == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap3 : ", *((ds1->getRMemoryPtr()->getVectorMemory())[0]) == *x0, true);
  delete ds1;
  cout << "--> swap test ended with success." << endl;
}


// plugins: plugins loading is already in testBuildDynamicalSystem2

void DynamicalSystemTest::End()
{
  cout << "==========================================" << endl;
  cout << " ===== End of DynamicalSystem tests ===== " << endl;
  cout << "==========================================" << endl;
}
