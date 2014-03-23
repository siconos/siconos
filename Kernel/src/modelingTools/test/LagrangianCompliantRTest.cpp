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
#include "LagrangianCompliantRTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianCompliantRTest);


void LagrangianCompliantRTest::setUp()
{
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("LagrangianCompliant_test.xml");
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

void LagrangianCompliantRTest::tearDown()
{}

// xml constructor (scleronomic case)
void LagrangianCompliantRTest::testBuildLagrangianCompliantR0()
{
  std::cout << "==============================================" <<std::endl;
  std::cout << "=== LagrangianCompliantR tests start ...=== " <<std::endl;
  std::cout << "==============================================" <<std::endl;
  SP::LagrangianCompliantR R1(new LagrangianCompliantR(tmpxml1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1a : ", R1->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1b : ", R1->getSubType() == RELATION::CompliantR, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1c : ", R1->gethName() == "TestPlugin:hCompl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1d : ", R1->getJachqName() == "TestPlugin:G0Compl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1e : ", R1->getJachlambdaName() == "TestPlugin:G1Compl", true);
  std::cout << " xml Constructor (1) LagrangianCompliantR ok" <<std::endl;
}

// data constructor:
void LagrangianCompliantRTest::testBuildLagrangianCompliantR2()
{
  SP::LagrangianCompliantR R1(new LagrangianCompliantR("TestPlugin:hCompl", "TestPlugin:G0Compl", "TestPlugin:G1Compl"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3a : ", R1->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3b : ", R1->getSubType() == RELATION::CompliantR, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3c : ", R1->gethName() == "TestPlugin:hCompl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3d : ", R1->getJachqName() == "TestPlugin:G0Compl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3e : ", R1->getJachlambdaName() == "TestPlugin:G1Compl", true);
  std::cout << " data Constructor LagrangianCompliantR ok" <<std::endl;
}


void LagrangianCompliantRTest::End()
{
  std::cout << "==============================================" <<std::endl;
  std::cout << " ===== End of LagrangianCompliantR tests ===== " <<std::endl;
  std::cout << "==============================================" <<std::endl;
}
