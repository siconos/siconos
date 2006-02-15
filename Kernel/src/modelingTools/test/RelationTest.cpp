/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
#include "RelationTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(RelationTest);


void RelationTest::setUp()
{
  // parse xml file:
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile("relation_test.xml");
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
  node1 = SiconosDOMTreeTools::findNodeChild(nodetmp, "Relation");
  tmpxml1 = new RelationXML(node1);
}

void RelationTest::tearDown()
{
  delete nsds;
  delete tmpxml1;
}

// default constructor
void RelationTest::testBuildRelation1()
{
  cout << "================================" << endl;
  cout << "=== Relation tests start ...=== " << endl;
  cout << "================================" << endl;
  Relation * R1 = new Relation();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation1a : ", R1->getInteractionPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation1b : ", R1->getType() == "Relation", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation1c : ", R1->getComputeOutputName() == "DefaultPlugin:computeOutput", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation1d : ", R1->getComputeInputName() == "DefaultPlugin:computeInput", true);
  delete R1;
  cout << " Default Constructor Relation ok" << endl;
}

// xml constructor
void RelationTest::testBuildRelation2()
{
  Interaction * inter = (nsds->getInteractions())[0];
  Relation * R1 = new Relation(tmpxml1, inter);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation2a : ", R1->getInteractionPtr()->getId() == "test-of-rel", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation2b : ", R1->getType() == "Relation", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation2c : ", R1->getComputeOutputName() == "TestPlugin:y", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation2d : ", R1->getComputeInputName() == "TestPlugin:R", true);
  delete R1;
  cout << " xml Constructor relation ok" << endl;
}
// copy constructor
void RelationTest::testBuildRelation3()
{
  Relation * R1 = new Relation();
  R1->setComputeInputFunction("TestPlugin.so", "R");
  R1->setComputeOutputFunction("TestPlugin.so", "y");

  Relation * R2 = new Relation(*R1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation3a : ", R2->getInteractionPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation3b : ", R2->getType() == "Relation", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation3c : ", R2->getComputeOutputName() == "TestPlugin:y", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildRelation3d : ", R2->getComputeInputName() == "TestPlugin:R", true);
  delete R2;
  delete R1;
  cout << " copy Constructor relation ok" << endl;
}


// computeOutput
/*void RelationTest::testComputeOutput()
{}
*/
// computeInput
/*void RelationTest::testComputeInput()
{}
*/
void RelationTest::End()
{
  cout << "====================================" << endl;
  cout << " ===== End of Relation Tests ===== " << endl;
  cout << "====================================" << endl;
}
