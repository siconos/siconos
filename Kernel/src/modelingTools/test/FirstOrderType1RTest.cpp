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

// \Warning these tests are not complete: add xml constructor.

#include "Interaction.h"
#include "FirstOrderType1RTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderType1RTest);


void FirstOrderType1RTest::setUp()
{
}

void FirstOrderType1RTest::tearDown()
{
}

// data constructor
void FirstOrderType1RTest::testBuildFirstOrderType1R1()
{
  cout << "======================================" << endl;
  cout << "=== FirstOrderType1R tests start ...== " << endl;
  cout << "======================================" << endl;
  cout << "--> Test: constructor data " << endl;

  FirstOrderType1R * R1 = new FirstOrderType1R("TestPlugin:hT1", "TestPlugin:gT1");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1a : ", R1->getInteractionPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1b : ", R1->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1c : ", R1->getSubType() == "Type1R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1d : ", R1->getFunctionName("h") == "TestPlugin:hT1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1e : ", R1->getFunctionName("g") == "TestPlugin:gT1", true);
  R1->setComputeJacobianHFunction("TestPlugin.so", "Jh0T1");
  R1->setComputeJacobianGFunction("TestPlugin.so", "Jg0T1");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1e : ", R1->getFunctionName("jacobianH0") == "TestPlugin:Jh0T1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1g : ", R1->getFunctionName("jacobianG0") == "TestPlugin:Jg0T1", true);
  delete R1;

  cout << "--> Constructor1 test ended with success." << endl;
}

void FirstOrderType1RTest::testBuildFirstOrderType1R2()
{
  cout << "--> Test: constructor data (2)." << endl;
  FirstOrderType1R * R2 = new FirstOrderType1R("TestPlugin:hT1", "TestPlugin:gT1", "TestPlugin:Jh0T1", "TestPlugin:Jg0T1");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2a : ", R2->getInteractionPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2b : ", R2->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2c : ", R2->getSubType() == "Type1R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2d : ", R2->getFunctionName("h") == "TestPlugin:hT1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2e : ", R2->getFunctionName("g") == "TestPlugin:gT1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2e : ", R2->getFunctionName("jacobianH0") == "TestPlugin:Jh0T1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2g : ", R2->getFunctionName("jacobianG0") == "TestPlugin:Jg0T1", true);
  delete R2;
  cout << "--> Constructor2 test ended with success." << endl;
}

// xml constructor
void FirstOrderType1RTest::testBuildFirstOrderType1R3()
{
  cout << "--> Test: constructor xml ." << endl;
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("FOT1R.xml");
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

  // look for NSDS, Interaction and relation node
  xmlNode* nodetmp = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  NonSmoothDynamicalSystemXML * nsdsxml = new NonSmoothDynamicalSystemXML(nodetmp);
  nsds = new NonSmoothDynamicalSystem(nsdsxml);
  delete nsdsxml;

  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Definition");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Content");

  // get relation
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderRelation");
  tmpxml1 = new FirstOrderRXML(node1);
  FirstOrderType1R * R1 = new FirstOrderType1R(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3a : ", R1->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3b : ", R1->getSubType() == "Type1R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3c : ", R1->getFunctionName("h") == "TestPlugin:hT1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3d : ", R1->getFunctionName("g") == "TestPlugin:gT1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3e : ", R1->getFunctionName("jacobianH0") == "TestPlugin:Jh0T1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3g : ", R1->getFunctionName("jacobianG0") == "TestPlugin:Jg0T1", true);
  delete R1;
  delete nsds;
  delete tmpxml1;
  cout << "--> Constructor xml test ended with success." << endl;
}

void FirstOrderType1RTest::End()
{
  cout << "==========================================" << endl;
  cout << " ===== End of FirstOrderType1R Tests ===== " << endl;
  cout << "==========================================" << endl;
}
