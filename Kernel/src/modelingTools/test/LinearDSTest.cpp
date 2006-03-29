/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#include "LinearDSTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LinearDSTest);


void LinearDSTest::setUp()
{
  x0 = new SimpleVector(3);
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  b0 = new SimpleVector(3);
  (*b0)(0) = 4;
  (*b0)(1) = 5;
  (*b0)(2) = 6;


  A0 = new SimpleMatrix("matA0.dat", true);

  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("lds_test.xml");
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LinearDS");
  tmpxml1 = new LinearDSXML(node1, false);
  // get second ds
  node2 = SiconosDOMTreeTools::findFollowNode(node1, "LinearDS");
  tmpxml2 = new LinearDSXML(node2, false);
}

void LinearDSTest::tearDown()
{
  delete tmpxml1;
  delete tmpxml2;
  delete x0;
  delete b0;
  delete A0;
}

// xml constructor (1), without plugin
void LinearDSTest::testBuildLinearDS1()
{
  cout << "======================================" << endl;
  cout << "=== LinearDS tests start ...=== " << endl;
  cout << "======================================" << endl;
  cout << "--> Test: constructor xml." << endl;
  LinearDS * ds = new LinearDS(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1A : ", ds->getType() == LDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1C : ", ds->getId() == "testDS1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1D : ", ds->getStepsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1E : ", ds->getX0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1F : ", ds->getB() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1G : ", ds->getA() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1H : ", ds->getIsLDSPlugin()[0] == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS1I : ", ds->getIsLDSPlugin()[1] == false, true);
  delete ds;
  cout << "--> Constructor xml test ended with success." << endl;
}

// xml constructor (2), with plugins
void LinearDSTest::testBuildLinearDS2()
{
  cout << "--> Test: constructor xml 2." << endl;
  LinearDS * ds = new LinearDS(tmpxml2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2A : ", ds->getType() == LDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2B : ", ds->getNumber() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2C : ", ds->getId() == "testDS2", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2E : ", ds->getX0() == 2 * *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2F : ", ds->getIsLDSPlugin()[0] == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2G : ", ds->getIsLDSPlugin()[1] == true, true);

  double time = 1.5;
  ds->computeB(time);
  ds->computeA(time);
  SimpleVector * x01 = new SimpleVector(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2I : ", ds->getB() == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2J : ", ds->getA() == 2 * *A0, true);

  SimpleVector* u0 = new SimpleVector(2);
  (*u0)(0) = 3;
  (*u0)(1) = 6;
  SiconosMatrix *T0 = new SimpleMatrix("matT0.dat", true);

  ds->computeVectorField(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2K : ", ds->getU() == *u0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2L : ", ds->getT() == *T0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2M : ", (*(ds->getXDotPtr())) == (*(ds->getAPtr()) * 2 * *x0 + ds->getB() + * (ds->getTPtr())**(ds->getUPtr())), true);

  delete T0;
  delete u0;
  delete x01;
  delete ds;
  cout << "--> Constructor xml 2 test ended with success." << endl;
}


// constructor from data
void LinearDSTest::testBuildLinearDS3()
{
  cout << "--> Test: constructor 3." << endl;
  LinearDS * ds = new LinearDS(13, 3, *x0, "TestPlugin:computeA", "TestPlugin:computeB");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3A : ", ds->getType() == LDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3C : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3D : ", ds->getX0() == *x0, true);

  double time = 1.5;
  ds->computeA(time);
  ds->computeB(time);
  ds->computeVectorField(time);
  SimpleVector * x01 = new SimpleVector(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3I : ", ds->getB() == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3J : ", ds->getA() == 2 * *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3E : ", ds->getXDot() == (time* *x01 + 2 * *A0 **x0) , true);

  ds->setUSize(3);
  ds->setComputeUFunction("TestPlugin.so", "computeU");
  ds->computeU(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3K : ", ds->getU() == time **x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS3L : ", ds->getTPtr() == NULL, true);
  ds->computeVectorField(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS2M : ", *(ds->getXDotPtr()) == (2 * *A0 * *x0 + ds->getB() + time**x0), true);
  delete ds;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// copy constructor
void LinearDSTest::testBuildLinearDS4()
{
  cout << "--> Test: constructor 4." << endl;
  LinearDS * ds1 = new LinearDS(tmpxml1);
  LinearDS * ds2 = new LinearDS(*ds1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4A : ", ds2->getType() == LDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4B : ", ds2->getNumber() == -2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4C : ", ds2->getId() == "copy", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4D : ", ds2->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4E : ", ds2->getX0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4F : ", ds2->getB() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4G : ", ds2->getA() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4H : ", ds2->getIsLDSPlugin()[0] == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearDS4I : ", ds2->getIsLDSPlugin()[1] == false, true);

  delete ds1;
  delete ds2;
  cout << "--> Constructor 4 test ended with success." << endl;
}

// setA
void LinearDSTest::testSetA()
{
  cout << "--> Test: setA." << endl;
  LinearDS * ds1 = new LinearDS(tmpxml2);
  ds1->setA(*A0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAa : ", ds1->getA() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAb : ", *(ds1->getJacobianXPtr()) == *A0, true);
  delete ds1;
  cout << "--> setA test ended with success." << endl;
}

// setAPtr
void LinearDSTest::testSetAPtr()
{
  cout << "--> Test: setAPtr." << endl;
  LinearDS * ds1 = new LinearDS(tmpxml2);
  ds1->setAPtr(A0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr : ", ds1->getAPtr() == A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAb : ", ds1->getJacobianXPtr() == A0, true);
  delete ds1;
  cout << "--> setAPtr test ended with success." << endl;
}


// setB
void LinearDSTest::testSetB()
{
  cout << "--> Test: setB." << endl;
  LinearDS * ds1 = new LinearDS(tmpxml2);
  ds1->setB(*b0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB : ", ds1->getB() == *b0, true);
  delete ds1;
  cout << "--> setB test ended with success." << endl;
}

// setBPtr
void LinearDSTest::testSetBPtr()
{
  cout << "--> Test: setBPtr." << endl;
  LinearDS* ds1 = new LinearDS(tmpxml2);
  ds1->setBPtr(b0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ds1->getBPtr() == b0, true);
  delete ds1;
  cout << "--> setBPtr test ended with success." << endl;
}

// plugins: plugins loading is already in testBuildLinearDS2

void LinearDSTest::End()
{
  cout << "==================================" << endl;
  cout << " ===== End of LinearDS tests ===== " << endl;
  cout << "==================================" << endl;
}
