/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#include "LagrangianDSTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianDSTest);


void LagrangianDSTest::setUp()
{
  q0 = new SimpleVector(3);
  (*q0)(0) = 1;
  (*q0)(1) = 2;
  (*q0)(2) = 3;
  velocity0 = new SimpleVector(3);
  (*velocity0)(0) = 4;
  (*velocity0)(1) = 5;
  (*velocity0)(2) = 6;

  u0 = new SimpleVector(2);
  (*u0)(0) = 4;
  (*u0)(1) = 5;

  mass = new SimpleMatrix(3, 3);
  (*mass)(0, 0) = 1;
  (*mass)(1, 1) = 2;
  (*mass)(2, 2) = 3;

  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("lagds_test.xml");
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LagrangianDS");
  tmpxml1 = new LagrangianDSXML(node1, false);
  // get second ds
  node2 = SiconosDOMTreeTools::findFollowNode(node1, "LagrangianDS");
  tmpxml2 = new LagrangianDSXML(node2, false);
  // get third ds
  node3 = SiconosDOMTreeTools::findFollowNode(node2, "LagrangianDS");
  tmpxml3 = new LagrangianDSXML(node3, false);
}

void LagrangianDSTest::tearDown()
{
  delete tmpxml1;
  delete tmpxml2;
  delete tmpxml3;
  delete q0;
  delete velocity0;
  delete u0;
  delete mass;
}

// xml constructor (1), without plugin
void LagrangianDSTest::testBuildLagrangianDS1()
{
  cout << "===================================" << endl;
  cout << "=== LagrangianDS tests start ...=== " << endl;
  cout << "===================================" << endl;
  cout << "--> Test: constructor xml." << endl;
  LagrangianDS * ds = new LagrangianDS(tmpxml1);
  SimpleMatrix M(3, 3);
  M.eye();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1A : ", ds->getType() == LNLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1C : ", ds->getId() == "testLAGDS1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1D : ", ds->getStepsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1E : ", ds->getQ0() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1F : ", ds->getVelocity0() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1E : ", ds->getQ() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1F : ", ds->getVelocity() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1F : ", ds->getMass() == M, true);
  delete ds;
  cout << "--> Constructor xml test ended with success." << endl;
}


// xml constructor (2), with plugins
void LagrangianDSTest::testBuildLagrangianDS2()
{
  cout << "--> Test: constructor xml 2." << endl;
  LagrangianDS * ds = new LagrangianDS(tmpxml2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2A : ", ds->getType() == LNLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2B : ", ds->getNumber() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2C : ", ds->getId() == "testLAGDS2", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2D : ", ds->getNdof() == 3, true);

  double time = 1.5;
  ds->initialize("TimeStepping", time);

  SimpleVector * x01 = new SimpleVector(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  SimpleVector * x02 = new SimpleVector(3);
  (*x02)(0) = 0;
  (*x02)(1) = 1 * (*q0)(1);
  (*x02)(2) = 2 * (*q0)(2);

  SimpleMatrix M(3, 3);
  M(0, 0) = 1;
  M(1, 1) = 2;
  M(2, 2) = 3;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2I : ", ds->getFExt() == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2J : ", ds->getFInt() == *x02, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2K : ", ds->getNNL() == *x02, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2L : ", ds->getMass() == M, true);
  M(0, 0) = 0;
  M(0, 1) = 3;
  M(0, 2) = 6;
  M(1, 0) = 1;
  M(1, 1) = 4;
  M(1, 2) = 7;
  M(2, 0) = 2;
  M(2, 1) = 5;
  M(2, 2) = 8;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2M : ", ds->getJacobianFInt(0) == M, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2N : ", ds->getJacobianFInt(1) == M, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2O : ", ds->getJacobianNNL(0) == M, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2P : ", ds->getJacobianNNL(1) == M, true);

  delete x01;
  delete x02;
  delete ds;
  cout << "--> Constructor xml 2 test ended with success." << endl;
}

// xml constructor (3)
void LagrangianDSTest::testBuildLagrangianDS3()
{
  cout << "--> Test: constructor xml 3." << endl;
  LagrangianDS * ds = new LagrangianDS(tmpxml3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2A : ", ds->getType() == LNLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2B : ", ds->getNumber() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2C : ", ds->getId() == "testLAGDS3", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS2D : ", ds->getNdof() == 3, true);
  double time = 1.5;
  ds->initialize("TimeStepping", time);

  SimpleMatrix M(3, 3);
  M(0, 0) = 1;
  M(1, 1) = 2;
  M(2, 2) = 3;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3I : ", ds->getFExt() == 2 * *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3J : ", ds->getFInt() == 3 * *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3K : ", ds->getNNL() ==  4 * *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3L : ", ds->getMass() == M, true);
  M.zero();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3M : ", ds->getJacobianFInt(0) == M, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3N : ", ds->getJacobianFInt(1) == M, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3O : ", ds->getJacobianNNL(0) == M, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3P : ", ds->getJacobianNNL(1) == M, true);

  map<string, bool> isPl = ds->getIsPlugin();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3Q : ", isPl["fExt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3R : ", isPl["fIxt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3S : ", isPl["NNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3T : ", isPl["jacobianQFInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3T : ", isPl["jacobianVelocityFInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3T : ", isPl["jacobianQNNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS3T : ", isPl["jacobianVelocityNNL"], false);

  delete ds;
  cout << "--> Constructor xml 3 test ended with success." << endl;
}

// constructor from data
void LagrangianDSTest::testBuildLagrangianDS4()
{
  cout << "--> Test: constructor 4." << endl;

  LagrangianDS * ds = new LagrangianDS(13, *q0, *velocity0, (*mass));
  double time = 1.5;
  ds->initialize("TimeStepping", time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4A : ", ds->getType() == LNLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4L : ", ds->getMass() == (*mass), true);

  map<string, bool> isPl = ds->getIsPlugin();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4Q : ", isPl["fExt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4R : ", isPl["fInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4S : ", isPl["NNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianQFInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianVelocityFInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianQNNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianVelocityNNL"], false);


  delete ds;
  cout << "--> Constructor 4 test ended with success." << endl;
}

// constructor from data
void LagrangianDSTest::testBuildLagrangianDS5()
{
  cout << "--> Test: constructor 5." << endl;
  string plugin = "TestPlugin:computeMass";
  DynamicalSystem * ds = new LagrangianDS(13, *q0, *velocity0, plugin);
  double time = 1.5;
  ds->initialize("TimeStepping", time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5A : ", ds->getType() == LNLDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5B : ", ds->getNumber() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5D : ", static_cast<LagrangianDS*>(ds)->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5L : ", static_cast<LagrangianDS*>(ds)->getMass() == (*mass), true);

  map<string, bool> isPl = ds->getIsPlugin();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5M : ", isPl["mass"], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5Q : ", isPl["fExt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5R : ", isPl["fInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5S : ", isPl["NNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianQFInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianVelocityFInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianQNNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianVelocityNNL"], false);


  delete ds;
  cout << "--> Constructor 5 test ended with success." << endl;
}

void LagrangianDSTest::testcomputeDS()
{
  cout << "-->Test: computeDS." << endl;
  DynamicalSystem * ds = new LagrangianDS(tmpxml2);
  LagrangianDS * copy = static_cast<LagrangianDS*>(ds);
  double time = 1.5;
  ds->initialize("EventDriven", time);
  ds->computeRhs(time);
  cout << "-->Test: computeDS." << endl;
  ds->computeJacobianXRhs(time);
  cout << "-->Test: computeDS." << endl;
  SimpleMatrix M(3, 3);
  M(0, 0) = 1;
  M(1, 1) = 2;
  M(2, 2) = 3;
  SiconosMatrix * jx = ds->getJacobianXRhsPtr();
  SiconosVector * vf = ds->getRhsPtr();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSI : ", *(vf->getVectorPtr(0)) == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSJ : ", prod(M, *(vf->getVectorPtr(1))) == (copy->getFExt() - copy->getFInt() - copy->getNNL()) , true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->getBlockPtr(1, 0))) == (copy->getJacobianFL(0)) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->getBlockPtr(1, 1))) == (copy->getJacobianFL(1)) , true);

  delete ds;
  cout << "--> computeDS test ended with success." << endl;


}
void LagrangianDSTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of LagrangianDS tests ===== " << endl;
  cout << "======================================" << endl;
}
