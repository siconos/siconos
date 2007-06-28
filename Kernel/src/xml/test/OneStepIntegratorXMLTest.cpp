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
#include "OneStepIntegratorXMLTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(OneStepIntegratorXMLTest);




void OneStepIntegratorXMLTest::setUp()
{
  try
  {
    modelXML = new SiconosModelXML("OneStepI.xml");
    strategyXML = modelXML->getStrategyXML();
    oneStepIs = strategyXML->getOneStepIntegratorXML();
    cout << "%%% setup : size == " << oneStepIs.size() << endl;
  }
  catch (SiconosException e)
  {
    cout << "Error in OneStepIntegratorXMLTest : " << e.report() << endl;
    exit(0);
  }
}

void OneStepIntegratorXMLTest::tearDown()
{

}

void OneStepIntegratorXMLTest::testNbOneStepIntegrator()
{
  CPPUNIT_ASSERT_MESSAGE("testNbOneStepIntegrator ", oneStepIs.size() == 3);
  cout << "OneStepIntegratorXMLTest >>> testNbOneStepIntegrator ................ OK\n ";
}

void OneStepIntegratorXMLTest::testGetR()
{
  // integrator 0 is a Moreau and doesn't have a attribute r !
  //  CPPUNIT_ASSERT_MESSAGE("testGetR : OneStepIntegator 1", oneStepIs[0]->getR() == 55);
  CPPUNIT_ASSERT_MESSAGE("testGetR : OneStepIntegator 2", oneStepIs[1]->getR() == 0);
  CPPUNIT_ASSERT_MESSAGE("testGetR : OneStepIntegator 3", oneStepIs[2]->getR() == -55);
  cout << "OneStepIntegratorXMLTest >>> testGetR ............................... OK\n ";
}


void OneStepIntegratorXMLTest::testGetDSConcerned()
{

  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : OneStepIntegator 1 :size", oneStepIs[0]->getDSConcerned().size() == 1);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : OneStepIntegator 1", oneStepIs[0]->getDSConcerned().at(0) == 2);

  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : OneStepIntegator 2 :size", oneStepIs[1]->getDSConcerned().size() == 1);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : OneStepIntegator 2", oneStepIs[1]->getDSConcerned().at(0) == 2);

  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : OneStepIntegator 3 :size", oneStepIs[2]->getDSConcerned().size() == 1);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : OneStepIntegator 3", oneStepIs[2]->getDSConcerned().at(0) == 1);

  cout << "OneStepIntegratorXMLTest >>> testGetDSConcerned ..................... OK\n ";
}



void OneStepIntegratorXMLTest::testGetType()
{
  CPPUNIT_ASSERT_MESSAGE("testGetType : OneStepIntegator 1 : Moreau ", oneStepIs[0]->getType() == "Moreau");
  CPPUNIT_ASSERT_MESSAGE("testGetType : OneStepIntegator 2 : Adams", oneStepIs[1]->getType() == "Adams");
  CPPUNIT_ASSERT_MESSAGE("testGetType : OneStepIntegator 3 : Lsodar", oneStepIs[2]->getType() == "LSODAR");

  cout << "OneStepIntegratorXMLTest >>> testGetType ............................ OK\n ";
}

void OneStepIntegratorXMLTest::testAdamsXML()
{
  AdamsXML* adams = static_cast<AdamsXML*>(oneStepIs[1]);

  CPPUNIT_ASSERT_MESSAGE("testAdamsXML type ", adams->getType() == "Adams");
  CPPUNIT_ASSERT_MESSAGE("testAdamsXML R", adams->getR() == 0);

  cout << "OneStepIntegratorXMLTest >>> testAdamsXML ........................... OK\n ";
}

void OneStepIntegratorXMLTest::testMoreauXML()
{
  MoreauXML* moreau = static_cast<MoreauXML*>(oneStepIs[0]);

  CPPUNIT_ASSERT_MESSAGE("testMoreauXML type ", moreau->getType() == "Moreau");
  //CPPUNIT_ASSERT_MESSAGE("testMoreauXML R", moreau->getR() == 55);

  cout << "OneStepIntegratorXMLTest >>> testMoreauXML .......................... OK\n ";
}

void OneStepIntegratorXMLTest::testLsodarXML()
{
  LsodarXML* lsodar = static_cast<LsodarXML*>(oneStepIs[2]);

  CPPUNIT_ASSERT_MESSAGE("testLsodarXML type ", lsodar->getType() == "LSODAR");
  CPPUNIT_ASSERT_MESSAGE("testLsodarXML R", lsodar->getR() == -55);

  cout << "OneStepIntegratorXMLTest >>> testLsodarXML .......................... OK\n ";
}
