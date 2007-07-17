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
#include "NonSmoothDynamicalSystemTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NonSmoothDynamicalSystemTest);


void NonSmoothDynamicalSystemTest::setUp()
{
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("nsds_test.xml");
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
  node = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  tmpxml = new NonSmoothDynamicalSystemXML(node);
}

void NonSmoothDynamicalSystemTest::tearDown()
{
  delete tmpxml;
}

// xml constructor
void NonSmoothDynamicalSystemTest::testBuildNonSmoothDynamicalSystem1()
{
  cout << "====================================================" << endl;
  cout << " ===== NonSmoothDynamicalSystem tests start ...===== " << endl;
  cout << "====================================================" << endl;
  cout << "------- Xml Constructor test -------" << endl;
  NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(tmpxml);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemA : ", nsds->getDSVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemB : ", nsds->getDynamicalSystemPtr(0)->getNumber() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemC : ", nsds->getDynamicalSystemPtr(1)->getNumber() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemD : ", nsds->getInteractionVectorSize() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemE : ", nsds->getInteractionPtr(0)->getNumber() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemF : ", nsds->isBVP() == false, true);
  delete nsds;
  cout << " ------- Constructor xml NonSmoothDynamicalSystem ok -------" << endl;
}
// copy constructor
void NonSmoothDynamicalSystemTest::testBuildNonSmoothDynamicalSystem2()
{
  cout << "------- Copy Constructor test -------" << endl;
  NonSmoothDynamicalSystem * nsds1 = new NonSmoothDynamicalSystem(tmpxml);
  cout << " okokok0 " << endl;
  NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(*nsds1);
  cout << " okokok " << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2A : ", nsds->getDSVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2B : ", nsds->getDynamicalSystemPtr(0)->getNumber() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2C : ", nsds->getDynamicalSystemPtr(1)->getNumber() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2D : ", nsds->getInteractionVectorSize() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->getInteractionPtr(0)->getNumber() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2F : ", nsds->isBVP() == false, true);
  delete nsds;
  delete nsds1;
  cout << "------- Constructor copy NonSmoothDynamicalSystem ok -------" << endl;
}


// addDynamicalSystem
void NonSmoothDynamicalSystemTest::testaddDynamicalSystem()
{
  NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(tmpxml);
  xmlNode *node2 = SiconosDOMTreeTools::findNodeChild(node, "DS_Definition");
  xmlNode * node3 = SiconosDOMTreeTools::findNodeChild(node2, "LagrangianLinearTIDS");
  DynamicalSystemXML * tmpdsxml = new LagrangianLinearTIDSXML(node3, false);

  DynamicalSystem * ltids = new LagrangianLinearTIDS(tmpdsxml);
  ltids ->setNumber(23);

  nsds->addDynamicalSystem(ltids);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testaddDynamicalSystemA : ", nsds->getDSVectorSize() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testaddDynamicalSystemB: ", nsds->getDynamicalSystemPtr(0)->getNumber() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testaddDynamicalSystemC : ", nsds->getDynamicalSystemPtr(1)->getNumber() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testaddDynamicalSystemC : ", nsds->getDynamicalSystemPtr(2)->getNumber() == 23, true);
  delete ltids;
  delete nsds;
  delete tmpdsxml;

  cout << "------- test addDynamicalSystem ok -------" << endl;
}

// addInteraction
void NonSmoothDynamicalSystemTest::testaddInteraction()
{
  NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(tmpxml);
  xmlNode *node2 = SiconosDOMTreeTools::findNodeChild(node, "Interaction_Definition");
  xmlNode *node3 = SiconosDOMTreeTools::findNodeChild(node2, "Interaction");
  vector<int> tmp;
  tmp.resize(2, 1);
  InteractionXML* interxml = new InteractionXML(node3, tmp);
  Interaction * inter = new Interaction(interxml);
  nsds->addInteraction(inter);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2D : ", nsds->getInteractionVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->getInteractionPtr(0)->getNumber() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->getInteractionPtr(1)->getNumber() == 12, true);
  delete nsds;
  delete inter;
  delete interxml;
  cout << " ------- test addInteractiontest ok -------" << endl;
}

void NonSmoothDynamicalSystemTest::End()
{
  cout << "===================================================" << endl;
  cout << " ===== End of NonSmoothDynamicalSystem tests ===== " << endl;
  cout << "===================================================" << endl;
}
