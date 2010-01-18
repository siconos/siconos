/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
#include "LagrangianScleronomousRTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianScleronomousRTest);


void LagrangianScleronomousRTest::setUp()
{
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("LagrangianScleronomous_test.xml");
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
  xmlNode * nodetmp = SiconosDOMTreeTools::findNodeChild(node, "Interaction");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Content");
  // get relation
  xmlNode * node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LagrangianRelation");
  tmpxml1.reset(new RelationXML(node1));
}

void LagrangianScleronomousRTest::tearDown()
{}

// xml constructor (scleronomic case)
void LagrangianScleronomousRTest::testBuildLagrangianScleronomousR0()
{
  cout << "==============================================" << endl;
  cout << "=== LagrangianScleronomousR tests start ...=== " << endl;
  cout << "==============================================" << endl;
  SP::LagrangianScleronomousR R1(new LagrangianScleronomousR(tmpxml1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR1a : ", R1->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR1b : ", R1->getSubType() == RELATION::ScleronomousR, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR1c : ", R1->gethName() == "TestPlugin:hSclero", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR1d : ", R1->getJachqName() == "TestPlugin:G0Sclero", true);
  cout << " xml Constructor (1) LagrangianScleronomousR ok" << endl;
}

// data constructor:
void LagrangianScleronomousRTest::testBuildLagrangianScleronomousR2()
{
  SP::LagrangianScleronomousR R1(new LagrangianScleronomousR("TestPlugin:hSclero", "TestPlugin:G0Sclero"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR3a : ", R1->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR3b : ", R1->getSubType() == RELATION::ScleronomousR, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR3c : ", R1->gethName() == "TestPlugin:hSclero", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianScleronomousR3d : ", R1->getJachqName() == "TestPlugin:G0Sclero", true);
  cout << " data Constructor LagrangianScleronomousR ok" << endl;
}


void LagrangianScleronomousRTest::End()
{
  cout << "=================================================" << endl;
  cout << " ===== End of LagrangianScleronomousR tests ===== " << endl;
  cout << "=================================================" << endl;
}
