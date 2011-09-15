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
#include "FirstOrderLinearDSTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearDSTest);


void FirstOrderLinearDSTest::setUp()
{
  x0.reset(new SimpleVector(3));
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  b0.reset(new SimpleVector(3));
  (*b0)(0) = 4;
  (*b0)(1) = 5;
  (*b0)(2) = 6;

  A0.reset(new SimpleMatrix("matA0.dat", true));

  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("lds_test.xml");
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
  // look for NSDS node
  xmlNodePtr nodetmp = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "DS_Definition");
  // get first ds
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderLinearDS");
  tmpxml1.reset(new FirstOrderLinearDSXML(node1, false));
  // get second ds
  node2 = SiconosDOMTreeTools::findFollowNode(node1, "FirstOrderLinearDS");
  tmpxml2.reset(new FirstOrderLinearDSXML(node2, false));
}

void FirstOrderLinearDSTest::tearDown()
{}

// xml constructor (1), without plugin
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS1()
{
  cout << "======================================" << endl;
  cout << "=== FirstOrderLinearDS tests start ...=== " << endl;
  cout << "======================================" << endl;
  cout << "--> Test: constructor xml." << endl;
  SP::FirstOrderLinearDS ds(new FirstOrderLinearDS(tmpxml1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1A : ", Type::value(*ds) == Type::FirstOrderLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1B : ", ds->number() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1D : ", ds->getStepsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1E : ", ds->getX0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1F : ", *(ds->b()) == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1G : ", *(ds->A()) == *A0, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1H : ", ds->A()->isPlugged(),false);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1I : ", ds->b()->isPlugged(),false);
  cout << "--> Constructor xml test ended with success." << endl;
}

// xml constructor (2), with plugins
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS2()
{
  cout << "--> Test: constructor xml 2." << endl;
  SP::FirstOrderLinearDS ds(new FirstOrderLinearDS(tmpxml2));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2A : ", Type::value(*ds) == Type::FirstOrderLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2B : ", ds->number() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2E : ", ds->getX0() == 2 * *x0, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2F : ", ds->A()->isPlugged(), true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2G : ", ds->b()->isPlugged(), true);

  double time = 1.5;
  ds->initialize(0, 0, time);
  ds->computeb(time);
  ds->computeA(time);
  SP::SimpleVector x01(new SimpleVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  ds->b()->display();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2I : ", *(ds->b()) == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2J : ", *(ds->A()) == 2 * *A0, true);

  SP::SimpleVector u0(new SimpleVector(2));
  (*u0)(0) = 3;
  (*u0)(1) = 6;
  SP::SiconosMatrix T0(new SimpleMatrix("matT0.dat", true));

  ds->computeRhs(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2M : ", (*(ds->rhs())) == prod(*(ds->A()), 2 * *x0) + * (ds->b()) , true);

  cout << "--> Constructor xml 2 test ended with success." << endl;
}


// constructor from data
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS3()
{
  cout << "--> Test: constructor 3." << endl;
  SP::FirstOrderLinearDS ds(new FirstOrderLinearDS(*x0, "TestPlugin:computeA", "TestPlugin:computeb"));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3A : ", Type::value(*ds) == Type::FirstOrderLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3B : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3C : ", ds->getX0() == *x0, true);

  double time = 1.5;
  ds->initialize(0, 0, time);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeRhs(time);
  SP::SimpleVector x01(new SimpleVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3D : ", *(ds->b()) == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3E : ", *(ds->A()) == 2 * *A0, true);
  //  ds->rhs()->display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS3F : ", *(ds->rhs()) == (time* *x01 + 2 * prod(*A0, *x0)) , true);

  ds->computeRhs(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2M : ", *(ds->rhs()) == (2 * prod(*A0 , *x0) + * (ds->b())), true);
  cout << "--> Constructor 3 test ended with success." << endl;
}

// setA
// void FirstOrderLinearDSTest::testSetA()
// {
//   cout << "--> Test: setA." << endl;
//   SP::FirstOrderLinearDS ds1(new FirstOrderLinearDS(tmpxml2));
//   ds1->setA(*A0p);
//   //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAa : ", ds1->getA()==*A0, true);
//   cout << "--> setA test ended with success." << endl;
// }

// setAPtr
void FirstOrderLinearDSTest::testSetAPtr()
{
  cout << "--> Test: setAPtr." << endl;
  SP::FirstOrderLinearDS ds1(new FirstOrderLinearDS(tmpxml2));
  ds1->setAPtr(A0p);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr : ", ds1->A() == A0p, true);
  cout << "--> setAPtr test ended with success." << endl;
}


// setB
// void FirstOrderLinearDSTest::testSetB()
// {
//   cout << "--> Test: setB." << endl;
//   SP::FirstOrderLinearDS ds1(new FirstOrderLinearDS(tmpxml2));
//   ds1->setB(*b0);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB : ", ds1->getB()==*b0, true);
//   cout << "--> setB test ended with success." << endl;
// }

// setBPtr
void FirstOrderLinearDSTest::testSetBPtr()
{
  cout << "--> Test: setBPtr." << endl;
  FirstOrderLinearDS* ds1(new FirstOrderLinearDS(tmpxml2));
  ds1->setb(b0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ds1->b() == b0, true);
  cout << "--> setBPtr test ended with success." << endl;
}

// plugins: plugins loading is already in testBuildFirstOrderLinearDS2

void FirstOrderLinearDSTest::End()
{
  cout << "============================================" << endl;
  cout << " ===== End of FirstOrderLinearDS tests =====" << endl;
  cout << "============================================" << endl;
}
