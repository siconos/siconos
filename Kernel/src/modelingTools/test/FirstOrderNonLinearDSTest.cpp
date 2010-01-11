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
#include "FirstOrderNonLinearDSTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderNonLinearDSTest);


void FirstOrderNonLinearDSTest::setUp()
{
  x0.reset(new SimpleVector(3));
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  J0.reset(new SimpleMatrix("matJ0.dat", true));
  M.reset(new SimpleMatrix("matM.dat", true));

  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("ds_test.xml");
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
  xmlNode* nodetmp = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "DS_Definition");
  // get first ds
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderNonLinearDS");
  tmpxml1.reset(new FirstOrderNonLinearDSXML(node1, false));
  // get second ds
  node2 = SiconosDOMTreeTools::findFollowNode(node1, "FirstOrderNonLinearDS");
  tmpxml2.reset(new FirstOrderNonLinearDSXML(node2, false));
}

void FirstOrderNonLinearDSTest::tearDown()
{}

// xml constructor (1), without plugin
void FirstOrderNonLinearDSTest::testBuildFirstOrderNonLinearDS1()
{
  cout << "======================================" << endl;
  cout << "=== FirstOrderNonLinearDS tests start ...=== " << endl;
  cout << "======================================" << endl;
  cout << "--> Test: constructor xml." << endl;
  SP::FirstOrderNonLinearDS ds(new FirstOrderNonLinearDS(tmpxml1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1A : ", ds->getType() == DS::FONLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1D : ", ds->getStepsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1E : ", ds->getX0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1E : ", ds->getM() == *M, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1H : ", ds->F()->isPlugged(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1I : ", ds->jacobianXF()->isPlugged(), false);

  cout << "--> Constructor xml test ended with success." << endl;
}


// xml constructor (2), with plugins
void FirstOrderNonLinearDSTest::testBuildFirstOrderNonLinearDS2()
{
  cout << "--> Test: constructor xml 2." << endl;
  SP::FirstOrderNonLinearDS ds(new FirstOrderNonLinearDS(tmpxml2));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2A : ", ds->getType() == DS::FONLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2D : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2E : ", ds->getX0() == 2 * *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2F : ", ds->F()->isPlugged(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2G : ", ds->jacobianXF()->isPlugged(), true);
  ds->initialize("TimeStepping", 0.5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2H : ", ds->getF() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2I : ", ds->getJacobianXF() == *J0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2I : ", !ds->M(), false);
  cout << "--> Constructor xml 2 test ended with success." << endl;
}


// constructor from data
void FirstOrderNonLinearDSTest::testBuildFirstOrderNonLinearDS3()
{
  cout << "--> Test: constructor 3." << endl;
  SP::FirstOrderNonLinearDS ds(new FirstOrderNonLinearDS(*x0, "TestPlugin:computeF", "TestPlugin:computeJacobianXF"));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3A : ", ds->getType() == DS::FONLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3C : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3D : ", ds->getX0() == *x0, true);
  double time = 1.5;
  ds->initialize("TimeStepping", time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3E : ", ds->getRhs() == time* *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2H : ", ds->getF() == time* *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2I : ", ds->getJacobianXF() == *J0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2I : ", !ds->M(), false);
  cout << "--> Constructor 3 test ended with success." << endl;
}

// setX0
void FirstOrderNonLinearDSTest::testSetX0()
{
  cout << "--> Test: setX0." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX0(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0 : ", ds1->getX0() == *x0, true);
  cout << "--> setX0 test ended with success." << endl;
}

// setX0Ptr
void FirstOrderNonLinearDSTest::testSetX0Ptr()
{
  cout << "--> Test: setX0Ptr." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX0Ptr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0Ptr : ", ds1->getX0() == *x0, true);
  cout << "--> setX0Ptr test ended with success." << endl;
}

// setX
void FirstOrderNonLinearDSTest::testSetX()
{
  cout << "--> Test: setX." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX : ", ds1->getX() == *x0, true);
  cout << "--> setX test ended with success." << endl;
}

// setXPtr
void FirstOrderNonLinearDSTest::testSetXPtr()
{
  cout << "--> Test: setXPtr." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setXPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXPtr : ", ds1->getX() == *x0, true);
  cout << "--> setXPtr test ended with success." << endl;
}

// setR
void FirstOrderNonLinearDSTest::testSetR()
{
  cout << "--> Test: setR." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setR(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetR : ", ds1->getR() == *x0, true);
  cout << "--> setR test ended with success." << endl;
}

// setRPtr
void FirstOrderNonLinearDSTest::testSetRPtr()
{
  cout << "--> Test: setRPtr." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setRPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetRPtr : ", ds1->getR() == *x0, true);
  cout << "--> setRPtr test ended with success." << endl;
}

// set JacobianX
void FirstOrderNonLinearDSTest::testSetJacobianXF()
{
  cout << "--> Test: setJacobianXF." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setJacobianXF(*J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianX : ", ds1->getJacobianXF() == *J0, true);
  cout << "--> setJacobianXF test ended with success." << endl;
}

// setJacobianXPtr
void FirstOrderNonLinearDSTest::testSetJacobianXFPtr()
{
  cout << "--> Test: setJacobianXFPtr." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setJacobianXFPtr(J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianXFPtr : ", ds1->getJacobianXF() == *J0, true);
  cout << "--> setJacobianXFPtr test ended with success." << endl;
}

// init
void FirstOrderNonLinearDSTest::testInitMemory()
{
  cout << "--> Test: initMemory." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->initMemory(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem1 : ", ds1->xMemory()->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem3 : ", ds1->rMemory()->getMemorySize() == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem4 : ", ds1->xMemory()->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem6 : ", ds1->rMemory()->getNbVectorsInMemory() == 0, true);
  cout << "--> initMemory test ended with success." << endl;
}


// swap
void FirstOrderNonLinearDSTest::testSwap()
{
  cout << "--> Test: swap." << endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX(*x0);
  ds1->setR(*x0);
  ds1->initMemory(1);
  ds1->swapInMemory();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap1 : ", *((ds1->xMemory()->getVectorMemory())[0]) == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap3 : ", *((ds1->rMemory()->getVectorMemory())[0]) == *x0, true);
  cout << "--> swap test ended with success." << endl;
}


// plugins: plugins loading is already in testBuildFirstOrderNonLinearDS2

void FirstOrderNonLinearDSTest::End()
{
  cout << "===============================================" << endl;
  cout << " ===== End of FirstOrderNonLinearDS tests =====" << endl;
  cout << "===============================================" << endl;
}
