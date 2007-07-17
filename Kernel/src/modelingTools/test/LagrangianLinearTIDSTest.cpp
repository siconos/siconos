/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "LagrangianLinearTIDSTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearTIDSTest);


void LagrangianLinearTIDSTest::setUp()
{
  q0 = new SimpleVector(3);
  (*q0)(0) = 1;
  (*q0)(1) = 2;
  (*q0)(2) = 3;
  velocity0 = new SimpleVector(3);
  (*velocity0)(0) = 4;
  (*velocity0)(1) = 5;
  (*velocity0)(2) = 6;

  mass = new SimpleMatrix(3, 3);
  (*mass)(0, 0) = 1;
  (*mass)(1, 1) = 2;
  (*mass)(2, 2) = 3;

  K = new SimpleMatrix("K.dat", true);
  C = new SimpleMatrix("C.dat", true);


  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("lagtids_test.xml");
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LagrangianLinearTIDS");
  tmpxml1 = new LagrangianLinearTIDSXML(node1, false);
}

void LagrangianLinearTIDSTest::tearDown()
{
  delete tmpxml1;
  delete q0;
  delete velocity0;
  delete mass;
  delete K;
  delete C;
}

// xml constructor (1), without plugin
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS1()
{
  cout << "===========================================" << endl;
  cout << "=== LagrangianLinearTIDS tests start ...===" << endl;
  cout << "===========================================" << endl;
  cout << "--> Test: constructor xml." << endl;
  LagrangianLinearTIDS * ds = new LagrangianLinearTIDS(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1A : ", ds->getType() == LLTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1C : ", ds->getId() == "testLAGTIDS1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1D : ", ds->getStepsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1E : ", ds->getQ0() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1F : ", ds->getVelocity0() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1G : ", ds->getQ() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1H : ", ds->getVelocity() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1I : ", ds->getMass() == *mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1J : ", ds->getK() == *K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1K : ", ds->getC() == *C, true);
  delete ds;
  cout << "--> Constructor xml test ended with success." << endl;
}


// Mass, K, C
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS2()
{
  cout << "--> Test: constructor 2." << endl;
  LagrangianLinearTIDS * ds = new LagrangianLinearTIDS(8, *q0, *velocity0, *mass, *K, *C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2A : ", ds->getType() == LLTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2B : ", ds->getNumber() == 8, true);
  cout << "--> Test: constructor 2." << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getStepsInMemory() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2E : ", ds->getQ0() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2F : ", ds->getVelocity0() == *velocity0, true);
  cout << "--> Test: constructor 2." << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2G : ", ds->getQ() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2H : ", ds->getVelocity() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getMass() == *mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2J : ", ds->getK() == *K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2K : ", ds->getC() == *C, true);
  ds->setComputeFExtFunction("TestPlugin.so", "computeFExt");

  double time = 1.5;
  ds->initialize("TimeStepping", time);

  SimpleVector * x01 = new SimpleVector(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getFExt() == time* *x01, true);
  delete x01;
  delete ds;
  cout << "--> Constructor 2 test ended with success." << endl;
}

// only  Mass
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS3()
{
  cout << "--> Test: constructor 3." << endl;
  LagrangianLinearTIDS * ds = new LagrangianLinearTIDS(8, *q0, *velocity0, *mass);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2A : ", ds->getType() == LLTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2B : ", ds->getNumber() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getStepsInMemory() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2E : ", ds->getQ0() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2F : ", ds->getVelocity0() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2G : ", ds->getQ() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2H : ", ds->getVelocity() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getMass() == *mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2J : ", ds->getKPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2K : ", ds->getCPtr() == NULL, true);

  ds->setComputeFExtFunction("TestPlugin.so", "computeFExt");

  double time = 1.5;
  ds->initialize("TimeStepping", time);

  SimpleVector * x01 = new SimpleVector(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getFExt() == time* *x01, true);
  delete x01;
  delete ds;
  cout << "--> Constructor 3 test ended with success." << endl;
}

void LagrangianLinearTIDSTest::testcomputeDS()
{
  cout << "-->Test: computeDS." << endl;
  DynamicalSystem * ds = new LagrangianLinearTIDS(tmpxml1);
  LagrangianLinearTIDS * copy = static_cast<LagrangianLinearTIDS*>(ds);
  double time = 1.5;
  ds->initialize("EventDriven", time);
  SiconosMatrix * jx = ds->getJacobianXRhsPtr();
  SiconosVector * vf = ds->getRhsPtr();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSI : ", *(vf->getVectorPtr(0)) == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSJ : ", prod(*mass, *(vf->getVectorPtr(1))) == (copy->getFExt() - prod(*K, *(copy->getQPtr())) - prod(*C, *(copy->getVelocityPtr()))) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(*mass, *(jx->getBlockPtr(1, 0))) == (-1.0 * *K) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(*mass, *(jx->getBlockPtr(1, 1))) == (-1.0 * *C) , true);

  delete ds;
  cout << "--> computeDS test ended with success." << endl;


}
void LagrangianLinearTIDSTest::End()
{
  cout << "==============================================" << endl;
  cout << " ===== End of LagrangianLinearTIDS tests =====" << endl;
  cout << "==============================================" << endl;
}
