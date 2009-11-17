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
  node = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  tmpxml.reset(new NonSmoothDynamicalSystemXML(node));
}

void NonSmoothDynamicalSystemTest::tearDown()
{}

// xml constructor
void NonSmoothDynamicalSystemTest::testBuildNonSmoothDynamicalSystem1()
{
  cout << "====================================================" << endl;
  cout << " ===== NonSmoothDynamicalSystem tests start ...===== " << endl;
  cout << "====================================================" << endl;
  cout << "------- Xml Constructor test -------" << endl;
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(tmpxml));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemA : ", nsds->getDSVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemB : ", nsds->dynamicalSystem(0)->number() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemC : ", nsds->dynamicalSystem(1)->number() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemD : ", nsds->getInteractionVectorSize() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemE : ", nsds->interaction(0)->number() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystemF : ", nsds->isBVP() == false, true);
  cout << " ------- Constructor xml NonSmoothDynamicalSystem ok -------" << endl;
}
// copy constructor
void NonSmoothDynamicalSystemTest::testBuildNonSmoothDynamicalSystem2()
{
  cout << "------- Copy Constructor test -------" << endl;
  SP::NonSmoothDynamicalSystem  nsds1(new NonSmoothDynamicalSystem(tmpxml));
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(*nsds1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2A : ", nsds->getDSVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2B : ", nsds->dynamicalSystem(0)->number() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2C : ", nsds->dynamicalSystem(1)->number() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2D : ", nsds->getInteractionVectorSize() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->interaction(0)->number() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2F : ", nsds->isBVP() == false, true);
  cout << "------- Constructor copy NonSmoothDynamicalSystem ok -------" << endl;
}


// insertDynamicalSystem
void NonSmoothDynamicalSystemTest::testinsertDynamicalSystem()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(tmpxml));
  xmlNode *node2 = SiconosDOMTreeTools::findNodeChild(node, "DS_Definition");
  xmlNode * node3 = SiconosDOMTreeTools::findNodeChild(node2, "LagrangianLinearTIDS");
  SP::DynamicalSystemXML tmpdsxml(new LagrangianLinearTIDSXML(node3, false));

  SP::DynamicalSystem ltids(new LagrangianLinearTIDS(tmpdsxml));
  ltids ->setNumber(23);

  nsds->insertDynamicalSystem(ltids);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testinsertDynamicalSystemA : ", nsds->getDSVectorSize() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testinsertDynamicalSystemB: ", nsds->dynamicalSystem(0)->number() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testinsertDynamicalSystemC : ", nsds->dynamicalSystem(1)->number() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testinsertDynamicalSystemC : ", nsds->dynamicalSystem(2)->number() == 23, true);
  cout << "------- test insertDynamicalSystem ok -------" << endl;
}

// insertInteraction
void NonSmoothDynamicalSystemTest::testinsertInteraction()
{
  SP::NonSmoothDynamicalSystem  nsds(new NonSmoothDynamicalSystem(tmpxml));
  xmlNode *node2 = SiconosDOMTreeTools::findNodeChild(node, "Interaction_Definition");
  xmlNode *node3 = SiconosDOMTreeTools::findNodeChild(node2, "Interaction");
  vector<int> tmp;
  tmp.resize(2, 1);
  SP::InteractionXML interxml(new InteractionXML(node3, tmp));
  SP::Interaction inter(new Interaction(interxml));
  nsds->insertInteraction(inter);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2D : ", nsds->getInteractionVectorSize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->interaction(0)->number() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothDynamicalSystem2E : ", nsds->interaction(1)->number() == 12, true);
  cout << " ------- test insertInteractiontest ok -------" << endl;
}

void NonSmoothDynamicalSystemTest::End()
{
  cout << "===================================================" << endl;
  cout << " ===== End of NonSmoothDynamicalSystem tests ===== " << endl;
  cout << "===================================================" << endl;
}
