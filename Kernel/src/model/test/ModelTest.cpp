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
#include "ModelTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ModelTest);


void ModelTest::setUp()
{
  t0 = 0.1;
  T = 10.0;
}

void ModelTest::tearDown()
{}

void ModelTest::testBuildModel0()
{
  cout << "=============================" << endl;
  cout << "=== Model tests start ...=== " << endl;
  cout << "=============================" << endl;
  cout << "--> Test: constructor 0." << endl;
  SP::Model M(new Model(t0));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getT0() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getFinalT() == -1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getCurrentTime() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", !M->getSimulationPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", !M->getNonSmoothDynamicalSystemPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", !M->getSiconosModelXMLPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getTitle() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getAuthor() == "nobody", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getDescription() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getDate() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel0 : ", M->getXmlSchema() == "none", true);
  cout << "--> Constructor 0 test ended with success." << endl;
}

void ModelTest::testBuildModel1()
{
  cout << "--> Test: constructor 1." << endl;
  SP::Model M(new Model(t0, T));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getT0() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getFinalT() == T, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getCurrentTime() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getSimulationPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getNonSmoothDynamicalSystemPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getSiconosModelXMLPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getTitle() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getAuthor() == "nobody", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getDescription() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getDate() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getXmlSchema() == "none", true);
  cout << "--> Constructor 1 test ended with success." << endl;
}

void ModelTest::testBuildModel2()
{
  cout << "--> Test: constructor 2." << endl;
  SP::Model M(new Model(t0, T, "myModel", "SiconosTeam", "Description", "Today", "XMLschema"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getT0() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getFinalT() == T, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getCurrentTime() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getSimulationPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getNonSmoothDynamicalSystemPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getSiconosModelXMLPtr(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getTitle() == "myModel", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getAuthor() == "SiconosTeam", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getDescription() == "Description", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getDate() == "Today", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getXmlSchema() == "XMLschema", true);
  cout << "--> Constructor 2 test ended with success." << endl;
}

// xml constructor
void ModelTest::testBuildModel3()
{
  cout << "--> Test: constructor xml." << endl;
  std::string xmlFile = "ModelXml_test.xml" ;

  SP::Model M(new Model(xmlFile));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getT0() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getFinalT() == T, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getCurrentTime() == t0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getSimulationPtr(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getNonSmoothDynamicalSystemPtr(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", !M->getSiconosModelXMLPtr(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getTitle() == "tryMxml", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getAuthor() == "SiconosTeam", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getDescription() == "Description", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getDate() == "Today", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1 : ", M->getXmlSchema() == "../../../config/xmlschema/SiconosModelSchema-V1.2.xsd", true);
  cout << "--> Constructor xml test ended with success." << endl;
}

