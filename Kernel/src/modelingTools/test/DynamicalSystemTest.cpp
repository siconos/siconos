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

  T0 = new SiconosMatrix("matT0.dat", true);
  J0 = new SiconosMatrix("matJ0.dat", true);

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
  cout << " ------ Constructor xml (1) DynamicalSystem ok ------ " << endl;
}


// xml constructor (2), with plugins
void DynamicalSystemTest::testBuildDynamicalSystem2()
{
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
  cout << " ------ Constructor xml (2) DynamicalSystem ok ------ " << endl;
}


// constructor from data
void DynamicalSystemTest::testBuildDynamicalSystem3()
{
  DynamicalSystem * ds = new DynamicalSystem(13, 3, *x0, "TestPlugin:vectorField");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3A : ", ds->getType() == NLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3C : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3D : ", ds->getX0() == *x0, true);

  double time = 1.5;
  ds->computeVectorField(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildDynamicalSystem3E : ", ds->getXDot() == time* *x0, true);
  delete ds;
  cout << " ------ Constructor from data DynamicalSystem ok ------ " << endl;
}

// copy constructor
void DynamicalSystemTest::testBuildDynamicalSystem4()
{
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
  cout << " ------ Constructor copy DynamicalSystem ok ------ " << endl;
}

// setX0
void DynamicalSystemTest::testSetX0()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX0(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0 : ", ds1->getX0() == *x0, true);
  delete ds1;
  cout << " ------ SetX0 ok ------ " << endl;
}

// setX0Ptr
void DynamicalSystemTest::testSetX0Ptr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX0Ptr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0Ptr : ", ds1->getX0() == *x0, true);
  delete ds1;
  cout << " ------ SetX0Ptr ok ------ " << endl;
}

// setX0 with exception
void DynamicalSystemTest::testSetX02()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX0(*u0);
}

// setX
void DynamicalSystemTest::testSetX()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX : ", ds1->getX() == *x0, true);
  delete ds1;
  cout << " ------ SetX ok ------ " << endl;
}

// setXPtr
void DynamicalSystemTest::testSetXPtr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXPtr : ", ds1->getX() == *x0, true);
  delete ds1;
  cout << " ------ SetXPtr ok ------ " << endl;
}

// setX with exception
void DynamicalSystemTest::testSetX2()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setX(*u0);
  delete ds1;
}

// setXDot
void DynamicalSystemTest::testSetXDot()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXDot(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXDot : ", ds1->getXDot() == *x0, true);
  delete ds1;
  cout << " ------ SetXDot ok ------ " << endl;
}

// setXDotPtr
void DynamicalSystemTest::testSetXDotPtr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXDotPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXDotPtr : ", ds1->getXDot() == *x0, true);
  delete ds1;
  cout << " ------ SetXDotPtr ok ------ " << endl;
}

// setXDot with exception
void DynamicalSystemTest::testSetXDot2()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXDot(*u0);
  delete ds1;
}

// setXFree
void DynamicalSystemTest::testSetXFree()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXFree(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXFree : ", ds1->getXFree() == *x0, true);
  delete ds1;
  cout << " ------ SetXFree ok ------ " << endl;
}

// setXFreePtr
void DynamicalSystemTest::testSetXFreePtr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXFreePtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXFreePtr : ", ds1->getXFree() == *x0, true);
  delete ds1;
  cout << " ------ SetXFreePtr ok ------ " << endl;
}

// setXFree with exception
void DynamicalSystemTest::testSetXFree2()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setXFree(*u0);
  delete ds1;
}

// setR
void DynamicalSystemTest::testSetR()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setR(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetR : ", ds1->getR() == *x0, true);
  delete ds1;
  cout << " ------ SetR ok ------ " << endl;
}

// setRPtr
void DynamicalSystemTest::testSetRPtr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setRPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetRPtr : ", ds1->getR() == *x0, true);
  delete ds1;
  cout << " ------ SetRPtr ok ------ " << endl;
}

// setR with exception
void DynamicalSystemTest::testSetR2()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setR(*u0);
  delete ds1;
}

// set JacobianX
void DynamicalSystemTest::testSetJacobianX()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setJacobianX(*J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianX : ", ds1->getJacobianX() == *J0, true);
  delete ds1;
  cout << " ------ SetJacobianX ok ------ " << endl;
}

// setJacobianXPtr
void DynamicalSystemTest::testSetJacobianXPtr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setJacobianXPtr(J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianXPtr : ", ds1->getJacobianX() == *J0, true);
  delete ds1;
  cout << " ------ SetJacobianXPtr ok ------ " << endl;
}

// setJacobianX with exception
void DynamicalSystemTest::testSetJacobianX2()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setJacobianX(*T0);
  delete ds1;
}

// setU
void DynamicalSystemTest::testSetU()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setU(*u0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetU : ", ds1->getU() == *u0, true);
  delete ds1;
  cout << " ------ SetU ok ------ " << endl;
}

// setUPtr
void DynamicalSystemTest::testSetUPtr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setUPtr(u0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetUPtr : ", ds1->getU() == *u0, true);
  delete ds1;
  cout << " ------ SetUPtr ok ------ " << endl;
}

// setU with exception
void DynamicalSystemTest::testSetU2()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setU(*x0);
  delete ds1;
}

// set T
void DynamicalSystemTest::testSetT()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setT(*T0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetT : ", ds1->getT() == *T0, true);
  delete ds1;
  cout << " ------ SetT ok ------ " << endl;
}

// setTPtr
void DynamicalSystemTest::testSetTPtr()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setTPtr(T0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetTPtr : ", ds1->getT() == *T0, true);
  delete ds1;
  cout << " ------ SetTPtr ok ------ " << endl;
}

// setT with exception
void DynamicalSystemTest::testSetT2()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->setT(*J0);
  delete ds1;
}

// init
void DynamicalSystemTest::testInitMemory()
{
  DynamicalSystem * ds1 = new DynamicalSystem(tmpxml2);
  ds1->initMemory(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem1 : ", ds1->getXMemoryPtr()->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem2 : ", ds1->getXDotMemoryPtr()->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem3 : ", ds1->getRMemoryPtr()->getMemorySize() == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem4 : ", ds1->getXMemoryPtr()->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem5 : ", ds1->getXDotMemoryPtr()->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem6 : ", ds1->getRMemoryPtr()->getNbVectorsInMemory() == 0, true);
  delete ds1;
  cout << " ------ SetInitMemory ok ------ " << endl;

}


// swap
void DynamicalSystemTest::testSwap()
{
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
  cout << " ------ SetSwap ok ------ " << endl;
}


// plugins: plugins loading is already in testBuildDynamicalSystem2

void DynamicalSystemTest::End()
{
  cout << "==========================================" << endl;
  cout << " ===== End of DynamicalSystem tests ===== " << endl;
  cout << "==========================================" << endl;
}
