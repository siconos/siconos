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

// \Warning these tests are not complete: add xml constructor.

#include "Interaction.hpp"
#include "FirstOrderType1RTest.hpp"
#include "SSLH.hpp"



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
  std::cout << "======================================" <<std::endl;
  std::cout << "=== FirstOrderType1R tests start ...== " <<std::endl;
  std::cout << "======================================" <<std::endl;
  std::cout << "--> Test: constructor 1 " <<std::endl;

  SP::FirstOrderType1R R1(new FirstOrderType1R("TestPlugin:hT1", "TestPlugin:gT1"));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1b : ", R1->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1c : ", R1->getSubType() == RELATION::Type1R, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1d : ", R1->gethName()=="TestPlugin:hT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1e : ", R1->getgName()=="TestPlugin:gT1", true);
  std::string plugin = "TestPlugin" + SSLH::getSharedLibraryExtension();
  R1->setComputeJachxFunction(plugin, "Jh0T1");
  R1->setComputeJacglambdaFunction(plugin, "Jg0T1");
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1e : ", R1->getJachName(0)=="TestPlugin:Jh0T1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R1g : ", R1->getJacgName(0)=="TestPlugin:Jg0T1", true);
  std::cout << "--> Constructor1 test ended with success." <<std::endl;
}

void FirstOrderType1RTest::testBuildFirstOrderType1R2()
{
  std::cout << "--> Test: constructor data (2)." <<std::endl;
  SP::FirstOrderType1R R2(new FirstOrderType1R("TestPlugin:hT1", "TestPlugin:gT1", "TestPlugin:Jh0T1", "TestPlugin:Jg0T1"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2b : ", R2->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2c : ", R2->getSubType() == RELATION::Type1R, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2d : ", R2->gethName()=="TestPlugin:hT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2e : ", R2->getgName()=="TestPlugin:gT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2e : ", R2->getJachName(0)=="TestPlugin:Jh0T1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R2g : ", R2->getJacgName(0)=="TestPlugin:Jg0T1", true);
  std::cout << "--> Constructor2 test ended with success." <<std::endl;
}

// xml constructor
void FirstOrderType1RTest::testBuildFirstOrderType1R3()
{
  std::cout << "--> Test: constructor xml ." <<std::endl;
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("FOT1R.xml");
  if (!doc)
    XMLException::selfThrow("Document not parsed successfully");
  cur = xmlDocGetRootElement(doc);
  if (!cur)
  {
    XMLException::selfThrow("empty document");
    xmlFreeDoc(doc);
  }
  // get rootNode
  std::cout << "--> Constructor xml test ended with success." <<std::endl;

  if (xmlStrcmp(cur->name, (const xmlChar *) "SiconosModel"))
  {
    XMLException::selfThrow("document of the wrong type, root node !=SiconosModel");
    xmlFreeDoc(doc);
  }
  std::cout << "--> Constructor xml test ended with success." <<std::endl;

  // look for NSDS, Interaction and relation node
  xmlNode* nodetmp = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
  SP::NonSmoothDynamicalSystemXML nsdsxml(new NonSmoothDynamicalSystemXML(nodetmp));
  nsds.reset(new NonSmoothDynamicalSystem(nsdsxml));
  std::cout << "--> Constructor xml test ended with success." <<std::endl;

  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Definition");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction");
  nodetmp = SiconosDOMTreeTools::findNodeChild(nodetmp, "Interaction_Content");

  // get relation
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "FirstOrderRelation");
  tmpxml1.reset(new RelationXML(node1));
  SP::FirstOrderType1R R1(new FirstOrderType1R(tmpxml1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3a : ", R1->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3b : ", R1->getSubType() == RELATION::Type1R, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3c : ", R1->gethName()=="TestPlugin:hT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3d : ", R1->getgName()=="TestPlugin:gT1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3e : ", R1->getJachName(0)=="TestPlugin:Jh0T1", true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderType1R3g : ", R1->getJacgName(0)=="TestPlugin:Jg0T1", true);
  std::cout << "--> Constructor xml test ended with success." <<std::endl;
}

void FirstOrderType1RTest::End()
{
  std::cout << "==========================================" <<std::endl;
  std::cout << " ===== End of FirstOrderType1R Tests ===== " <<std::endl;
  std::cout << "==========================================" <<std::endl;
}
