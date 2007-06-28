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
#include "ModelTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ModelTest);


void ModelTest::setUp()
{
  t0 = 0.0;
  T = 10.0;
}

void ModelTest::tearDown()
{}

// xml constructor
void ModelTest::testBuildModel1()
{
  cout << "=============================" << endl;
  cout << "=== Model tests start ...=== " << endl;
  cout << "=============================" << endl;
  char * xmlFile = "ModelXml_test.xml" ;

  cout << "--> Test: constructor xml (1)." << endl;
  Model * M = new Model(xmlFile);
  //check that xml link is not NULL
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1a : ", M->getSiconosModelXMLPtr()!=NULL, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1b : ", M->getT0()==0.1, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1c : ", M->getFinalT()==10.0, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1d : ", M->getTitle()=="tryMxml", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1e : ", M->getAuthor()=="SiconosTeam", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1f : ", M->getDate()=="Today", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1g : ", M->getXmlSchema()=="/share/SICONOS/SiconosModelSchema-V1.2.xsd", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1h : ", M->getDescription()=="none", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1i : ", M->getNonSmoothDynamicalSystemPtr()!=NULL, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1j : ", M->getSimulationPtr()!=NULL, true);
  delete M;
  cout << "--> Constructor xml test ended with success." << endl;
}

// data constructor (1)
void ModelTest::testBuildModel2()
{
  cout << "--> Test: constructor xml (2)." << endl;
  Model * M = new Model(t0, T, "tryMxml", "SiconosTeam", "none", "Today", "/share/SICONOS/SiconosModelSchema-V1.2.xsd");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2a : ", M->getSiconosModelXMLPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2b : ", M->getT0() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2c : ", M->getFinalT() == 10.0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2d : ", M->getTitle() == "tryMxml", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2e : ", M->getAuthor() == "SiconosTeam", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2f : ", M->getDate() == "Today", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2g : ", M->getXmlSchema() == "/share/SICONOS/SiconosModelSchema-V1.2.xsd", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2h : ", M->getDescription() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2i : ", M->getNonSmoothDynamicalSystemPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2j : ", M->getSimulationPtr() == NULL, true);
  delete M;
  cout << "--> Constructor xml (2) test ended with success." << endl;
}

// setNonSmoothDynamicalSystem
void ModelTest::testsetNonSmoothDynamicalSystemPtr()
{
  // Change since NonSmoothDynamicalSystem default constructor is now private. To be updated.
  //   NonSmoothDynamicalSystem * newNsds = new NonSmoothDynamicalSystem();
  //   Model * M = new Model(t0,T);
  //   M->setNonSmoothDynamicalSystemPtr(newNsds);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("test setNonSmoothDynamicalSystemPtr : ", M->getNonSmoothDynamicalSystemPtr()==newNsds, true);
  //   delete M;
  //   cout << " test setNonSmoothDynamicalSystemPtr ok " << endl;
}

// setSiconosModelXML
void ModelTest::testsetSiconosModelXMLPtr()
{
  cout << "--> Test: setSiconosModelXMLPtr." << endl;
  SiconosModelXML * newMxml = new SiconosModelXML();
  Model * M = new Model(t0, T);
  M->setSiconosModelXMLPtr(newMxml);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" test setSiconosModelXMLPtr: ", M->getSiconosModelXMLPtr() == newMxml, true);
  delete M;
  delete newMxml;
  cout << "--> setSiconosModelXMLPtr test ended with success." << endl;
}
void ModelTest::End()
{
  cout << "================================" << endl;
  cout << " ===== End of Model Tests ===== " << endl;
  cout << "================================" << endl;
}
