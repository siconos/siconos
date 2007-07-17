/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "FirstOrderRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderRTest);


void FirstOrderRTest::setUp()
{
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("firstOrderR_test.xml");
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
}

void FirstOrderRTest::tearDown()
{
  delete nsds;
  delete tmpxml1;
}

// data constructor
void FirstOrderRTest::testBuildFirstOrderR1()
{
  cout << "===================================" << endl;
  cout << "=== FirstOrderR tests start ...=== " << endl;
  cout << "===================================" << endl;
  FirstOrderR * R1 = new FirstOrderR("TestPlugin:y", "TestPlugin:R");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1a : ", R1->getInteractionPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1b : ", R1->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1c : ", R1->getSubType() == "R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1d : ", R1->getFunctionName("h") == "TestPlugin:y", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1e : ", R1->getFunctionName("g") == "TestPlugin:R", true);
  R1->setComputeJacobianHFunction("TestPlugin.so", "Jh0", 0);
  R1->setComputeJacobianHFunction("TestPlugin.so", "Jh1", 1);
  R1->setComputeJacobianGFunction("TestPlugin.so", "Jg0");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1e : ", R1->getFunctionName("jacobianH0") == "TestPlugin:Jh0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1f : ", R1->getFunctionName("jacobianH1") == "TestPlugin:Jh1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR1g : ", R1->getFunctionName("jacobianG0") == "TestPlugin:Jg0", true);
  delete R1;

  cout << "--> Constructor1 test ended with success." << endl;
}

// xml constructor
void FirstOrderRTest::testBuildFirstOrderR2()
{
  cout << "--> Test: constructor xml ." << endl;
  FirstOrderR * R1 = new FirstOrderR(tmpxml1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR2a : ", R1->getType() == "FirstOrder", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR2b : ", R1->getSubType() == "R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR2c : ", R1->getFunctionName("h") == "TestPlugin:y", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR2d : ", R1->getFunctionName("g") == "TestPlugin:R", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR2e : ", R1->getFunctionName("jacobianH0") == "TestPlugin:Jh0", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR2f : ", R1->getFunctionName("jacobianH1") == "TestPlugin:Jh1", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderR2g : ", R1->getFunctionName("jacobianG0") == "TestPlugin:Jg0", true);
  delete R1;
  cout << "--> Constructor xml test ended with success." << endl;
}

void FirstOrderRTest::End()
{
  cout << "============================================" << endl;
  cout << " ===== End of FirstOrderR Tests ===== " << endl;
  cout << "=============================================" << endl;
}
