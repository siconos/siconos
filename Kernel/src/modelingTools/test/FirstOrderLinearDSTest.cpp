/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
#include "FirstOrderLinearDSTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearDSTest);


void FirstOrderLinearDSTest::setUp()
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderLinearDS");
  tmpxml1 = new FirstOrderLinearDSXML(node1, false);
  // get second ds
  node2 = SiconosDOMTreeTools::findFollowNode(node1, "FirstOrderLinearDS");
  tmpxml2 = new FirstOrderLinearDSXML(node2, false);
}

void FirstOrderLinearDSTest::tearDown()
{
  delete tmpxml1;
  delete tmpxml2;
  delete x0;
  delete b0;
  delete A0;
}

// xml constructor (1), without plugin
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS1()
{
  cout << "======================================" << endl;
  cout << "=== FirstOrderLinearDS tests start ...=== " << endl;
  cout << "======================================" << endl;
  cout << "--> Test: constructor xml." << endl;
  FirstOrderLinearDS * ds = new FirstOrderLinearDS(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1A : ", ds->getType() == FOLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1C : ", ds->getId() == "testDS1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1D : ", ds->getStepsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1E : ", ds->getX0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1F : ", ds->getB() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1G : ", ds->getA() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1H : ", ds->isPlugged("A"), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1I : ", ds->isPlugged("b"), false);
  delete ds;
  cout << "--> Constructor xml test ended with success." << endl;
}

// xml constructor (2), with plugins
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS2()
{
  cout << "--> Test: constructor xml 2." << endl;
  FirstOrderLinearDS * ds = new FirstOrderLinearDS(tmpxml2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2A : ", ds->getType() == FOLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2B : ", ds->getNumber() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2C : ", ds->getId() == "testDS2", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2E : ", ds->getX0() == 2 * *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2F : ", ds->isPlugged("A"), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2G : ", ds->isPlugged("b"), true);

  double time = 1.5;
  ds->initialize("TimeStepping", time);
  //   ds->computeB(time);
  //   ds->computeA(time);
  SimpleVector * x01 = new SimpleVector(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2I : ", ds->getB() == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2J : ", ds->getA() == 2 * *A0, true);

  SimpleVector* u0 = new SimpleVector(2);
  (*u0)(0) = 3;
  (*u0)(1) = 6;
  SiconosMatrix *T0 = new SimpleMatrix("matT0.dat", true);

  ds->computeRhs(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2M : ", (*(ds->getRhsPtr())) == prod(*(ds->getAPtr()), 2 * *x0) + ds->getB() , true);

  delete T0;
  delete u0;
  delete x01;
  delete ds;
  cout << "--> Constructor xml 2 test ended with success." << endl;
}


// constructor from data
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS3()
{
  cout << "--> Test: constructor 3." << endl;
  FirstOrderLinearDS * ds = new FirstOrderLinearDS(13, *x0, "TestPlugin:computeA", "TestPlugin:computeB");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3A : ", ds->getType() == FOLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3C : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3D : ", ds->getX0() == *x0, true);

  double time = 1.5;
  ds->initialize("TimeStepping", time);
  //   ds->computeA(time);
  //   ds->computeB(time);
  //   ds->computeRhs(time);
  SimpleVector * x01 = new SimpleVector(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3I : ", ds->getB() == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3J : ", ds->getA() == 2 * *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3E : ", ds->getRhs() == (time* *x01 + 2 * prod(*A0, *x0)) , true);

  ds->computeRhs(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2M : ", *(ds->getRhsPtr()) == (2 * prod(*A0 , *x0) + ds->getB()), true);
  delete ds;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// setA
void FirstOrderLinearDSTest::testSetA()
{
  cout << "--> Test: setA." << endl;
  FirstOrderLinearDS * ds1 = new FirstOrderLinearDS(tmpxml2);
  ds1->setA(*A0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAa : ", ds1->getA() == *A0, true);
  delete ds1;
  cout << "--> setA test ended with success." << endl;
}

// setAPtr
void FirstOrderLinearDSTest::testSetAPtr()
{
  cout << "--> Test: setAPtr." << endl;
  FirstOrderLinearDS * ds1 = new FirstOrderLinearDS(tmpxml2);
  ds1->setAPtr(A0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr : ", ds1->getAPtr() == A0, true);
  delete ds1;
  cout << "--> setAPtr test ended with success." << endl;
}


// setB
void FirstOrderLinearDSTest::testSetB()
{
  cout << "--> Test: setB." << endl;
  FirstOrderLinearDS * ds1 = new FirstOrderLinearDS(tmpxml2);
  ds1->setB(*b0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB : ", ds1->getB() == *b0, true);
  delete ds1;
  cout << "--> setB test ended with success." << endl;
}

// setBPtr
void FirstOrderLinearDSTest::testSetBPtr()
{
  cout << "--> Test: setBPtr." << endl;
  FirstOrderLinearDS* ds1 = new FirstOrderLinearDS(tmpxml2);
  ds1->setBPtr(b0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ds1->getBPtr() == b0, true);
  delete ds1;
  cout << "--> setBPtr test ended with success." << endl;
}

// plugins: plugins loading is already in testBuildFirstOrderLinearDS2

void FirstOrderLinearDSTest::End()
{
  cout << "============================================" << endl;
  cout << " ===== End of FirstOrderLinearDS tests =====" << endl;
  cout << "============================================" << endl;
}
