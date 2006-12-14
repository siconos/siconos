/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "LagrangianRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianRTest);


void LagrangianRTest::setUp()
{
  G0 = new SimpleMatrix("matG0.dat", true);
  G1 = new SimpleMatrix("matG.dat", true);
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("LagrangianR_test.xml");
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

  // second file
  // parse xml file:
  xmlDocPtr doc2;
  xmlNodePtr cur2;
  doc2 = xmlParseFile("LagrangianR_test2.xml");
  if (doc2 == NULL)
    XMLException::selfThrow("Document not parsed successfully");
  cur2 = xmlDocGetRootElement(doc2);
  if (cur2 == NULL)
  {
    XMLException::selfThrow("empty document");
    xmlFreeDoc(doc2);
  }

  // get rootNode

  if (xmlStrcmp(cur2->name, (const xmlChar *) "SiconosModel"))
  {
    XMLException::selfThrow("document of the wrong type, root node !=SiconosModel");
    xmlFreeDoc(doc2);
  }
  // look for NSDS node
  node = SiconosDOMTreeTools::findNodeChild(cur2, "NSDS");
  xmlNode * nodetmp2 = SiconosDOMTreeTools::findNodeChild(node, "Interaction");
  nodetmp2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "Interaction_Content");
  // get relation
  xmlNode * node2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "LagrangianRelation");
  tmpxml2 = new LagrangianRXML(node2);

  // third file
  // parse xml file:
  doc2 = xmlParseFile("LagrangianR_test3.xml");
  if (doc2 == NULL)
    XMLException::selfThrow("Document not parsed successfully");
  cur2 = xmlDocGetRootElement(doc2);
  if (cur2 == NULL)
  {
    XMLException::selfThrow("empty document");
    xmlFreeDoc(doc2);
  }

  // get rootNode

  if (xmlStrcmp(cur2->name, (const xmlChar *) "SiconosModel"))
  {
    XMLException::selfThrow("document of the wrong type, root node !=SiconosModel");
    xmlFreeDoc(doc2);
  }
  // look for NSDS node
  node = SiconosDOMTreeTools::findNodeChild(cur2, "NSDS");
  nodetmp2 = SiconosDOMTreeTools::findNodeChild(node, "Interaction");
  nodetmp2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "Interaction_Content");
  // get relation
  node2 = SiconosDOMTreeTools::findNodeChild(nodetmp2, "LagrangianRelation");
  tmpxml3 = new LagrangianRXML(node2);
}

void LagrangianRTest::tearDown()
{
  delete G0;
  delete G1;
  delete tmpxml1;
  delete tmpxml2;
  delete tmpxml3;
}

// xml constructor (scleronomic case)
void LagrangianRTest::testBuildLagrangianR0()
{
  cout << "==================================" << endl;
  cout << "=== LagrangianR tests start ...=== " << endl;
  cout << "==================================" << endl;
  LagrangianR * R1 = new LagrangianR(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1a : ", R1->getType() == "LagrangianR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1b : ", R1->getHFunctionName() == "TestPlugin:h0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1c : ", R1->getGFunctionName() == "TestPlugin:G0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1a : ", R1->getLagrangianRelationType() == "scleronomic", true);
  delete R1;
  cout << " xml Constructor (1) LagrangianR ok" << endl;
}

// xml constructor (rhenomorous case)
void LagrangianRTest::testBuildLagrangianR1()
{
  LagrangianR * R1 = new LagrangianR(tmpxml2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1a : ", R1->getType() == "LagrangianR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1b : ", R1->getHFunctionName() == "TestPlugin:h1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1c : ", R1->getGFunctionName(0) == "TestPlugin:G10", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1c : ", R1->getGFunctionName(1) == "TestPlugin:G11", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1a : ", R1->getLagrangianRelationType() == "rhenomorous", true);
  delete R1;
  cout << " xml Constructor (2) LagrangianR ok" << endl;
}

// xml constructor (scleronomic+lambda case)
void LagrangianRTest::testBuildLagrangianR4()
{
  LagrangianR * R1 = new LagrangianR(tmpxml3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR4a : ", R1->getType() == "LagrangianR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR4b : ", R1->getHFunctionName() == "TestPlugin:h2", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR4c : ", R1->getGFunctionName(0) == "TestPlugin:G20", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR4d : ", R1->getGFunctionName(1) == "TestPlugin:G21", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR4e : ", R1->getLagrangianRelationType() == "scleronomic+lambda", true);
  delete R1;
  cout << " xml Constructor (2) LagrangianR ok" << endl;
}

// data constructor:
void LagrangianRTest::testBuildLagrangianR2()
{
  vector<string> G;
  G.reserve(1);
  G.push_back("TestPlugin:G0");
  LagrangianR * R1 = new LagrangianR("scleronomic", "TestPlugin:h0", G);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR2a : ", R1->getType() == "LagrangianR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR2b : ", R1->getHFunctionName() == "TestPlugin:h0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR2c : ", R1->getGFunctionName() == "TestPlugin:G0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR2a : ", R1->getLagrangianRelationType() == "scleronomic", true);
  delete R1;
  cout << " data Constructor LagrangianR ok" << endl;
}

// copy constructor
void LagrangianRTest::testBuildLagrangianR3()
{

  Relation * ref = new LagrangianR(tmpxml1);
  LagrangianR * R1 = new LagrangianR(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1a : ", R1->getType() == "LagrangianR", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1b : ", R1->getHFunctionName() == "TestPlugin:h0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1c : ", R1->getGFunctionName() == "TestPlugin:G0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianR1a : ", R1->getLagrangianRelationType() == "scleronomic", true);
  delete R1;
  delete ref;
  cout << "copy Constructor LagrangianR ok" << endl;
}


void LagrangianRTest::End()
{
  cout << "============================================" << endl;
  cout << " ===== End of LagrangianR tests ===== " << endl;
  cout << "============================================" << endl;
}
