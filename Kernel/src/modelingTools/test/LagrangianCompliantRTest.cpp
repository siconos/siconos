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
#include "LagrangianCompliantRTest.h"
using namespace std;

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
  xmlNode * nodetmp = SiconosDOMTreeTools::findNodeChild(node, "Interaction");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Content");
  // get relation
  xmlNode * node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "LagrangianRelation");
  tmpxml1 = new LagrangianRXML(node1);
}

void LagrangianCompliantRTest::tearDown()
{
  delete tmpxml1;
}

// xml constructor (scleronomic case)
void LagrangianCompliantRTest::testBuildLagrangianCompliantR0()
{
  cout << "==============================================" << endl;
  cout << "=== LagrangianCompliantR tests start ...=== " << endl;
  cout << "==============================================" << endl;
  LagrangianCompliantR * R1 = new LagrangianCompliantR(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1a : ", R1->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1b : ", R1->getSubType() == "CompliantR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1c : ", R1->getFunctionName("h") == "TestPlugin:hCompl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1d : ", R1->getFunctionName("G0") == "TestPlugin:G0Compl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR1e : ", R1->getFunctionName("G1") == "TestPlugin:G1Compl", true);
  delete R1;
  cout << " xml Constructor (1) LagrangianCompliantR ok" << endl;
}

// default constructor:
void LagrangianCompliantRTest::testBuildLagrangianCompliantR1()
{
  LagrangianCompliantR * R1 = new LagrangianCompliantR();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR2a : ", R1->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR2b : ", R1->getSubType() == "CompliantR", true);
  delete R1;
  cout << " default Constructor LagrangianCompliantR ok" << endl;
}

// data constructor:
void LagrangianCompliantRTest::testBuildLagrangianCompliantR2()
{
  vector<string> G;
  G.reserve(1);
  G.push_back("TestPlugin:G0Compl");
  G.push_back("TestPlugin:G1Compl");

  LagrangianCompliantR * R1 = new LagrangianCompliantR("TestPlugin:hCompl", G);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3a : ", R1->getType() == "Lagrangian", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3b : ", R1->getSubType() == "CompliantR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3c : ", R1->getFunctionName("h") == "TestPlugin:hCompl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3d : ", R1->getFunctionName("G0") == "TestPlugin:G0Compl", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantR3d : ", R1->getFunctionName("G1") == "TestPlugin:G1Compl", true);
  delete R1;
  cout << " data Constructor LagrangianCompliantR ok" << endl;
}


void LagrangianCompliantRTest::End()
{
  cout << "==============================================" << endl;
  cout << " ===== End of LagrangianCompliantR tests ===== " << endl;
  cout << "==============================================" << endl;
}
