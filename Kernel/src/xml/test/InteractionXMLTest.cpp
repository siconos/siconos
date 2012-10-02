/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "InteractionXMLTest.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(InteractionXMLTest);




void InteractionXMLTest::setUp()
{
  try
  {
    this->doc = xmlParseFile("Interaction.xml");
    this->root = xmlDocGetRootElement(doc);
    vector<int> DSConcerned(3);
    DSConcerned.at(0) = 1;
    DSConcerned.at(1) = 2;
    DSConcerned.at(2) = 3;

    this->interaction = InteractionXML(root, DSConcerned);

  }
  catch (SiconosException e)
  {
    cout << "Error in setup of InteractionXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void InteractionXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
  //doc = NULL;
  //root = NULL;
  //child = NULL;
}

void InteractionXMLTest::testGetNumber()
{
  CPPUNIT_ASSERT_MESSAGE("testGetNumber ", interaction.number() == 120);
  cout << "InteractionXMLTest >>> testGetNumber ................................ OK\n ";
}

void InteractionXMLTest::testGetId()
{
  CPPUNIT_ASSERT_MESSAGE("testGetId ", interaction.getId() == "Interaction1");
  cout << "InteractionXMLTest >>> testGetId .................................... OK\n ";
}

void InteractionXMLTest::testHasYLambda()
{
  CPPUNIT_ASSERT_MESSAGE("testHasYLambda : HasY ", interaction.hasY() == false);
  CPPUNIT_ASSERT_MESSAGE("testHasYLambda : HasLambda ", interaction.hasLambda() == false);
  cout << "InteractionXMLTest >>> testHasYLambda ............................... OK\n ";
}

void InteractionXMLTest::testGetDSConcerned()
{
  vector < vector <int> > v;
  v = interaction.getDSConcerned();
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : v.size() ", v.size() == 2);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : v.at(0).size() ", v.at(0).size() == 2);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : v.at(1).size() ", v.at(1).size() == 2);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : v.at(0).at(0) ", v.at(0).at(0) == 1);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : v.at(0).at(1) ", v.at(0).at(1) == 2);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : v.at(1).at(0) ", v.at(1).at(0) == 1);
  CPPUNIT_ASSERT_MESSAGE("testGetDSConcerned : v.at(1).at(1) ", v.at(1).at(1) == 3);
  cout << "InteractionXMLTest >>> testGetDSConcerned ........................... OK\n ";
}






