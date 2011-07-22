/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
//$id$

#include "SiconosMemoryXMLTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryXMLTest);


void SiconosMemoryXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("SiconosMemory.xml");
    this->root = xmlDocGetRootElement(doc);
    //    this->smxml = SiconosMemoryXML(root);
  }
  catch (SiconosException e)
  {
    cout << "Error in SiconosMemoryXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void SiconosMemoryXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void SiconosMemoryXMLTest::testHasMemory()
{
  this->smxml = new SiconosMemoryXML(SiconosDOMTreeTools::findNodeChild(this->root, "test1"));
  CPPUNIT_ASSERT_MESSAGE("test hasMemory == true", this->smxml->hasMemory() == true);
  delete smxml;
  cout << "SiconosMemoryXMLTest::testHasMemory - test1 .........................OK" << endl;


  this->smxml = new SiconosMemoryXML(SiconosDOMTreeTools::findNodeChild(this->root, "test2"));
  CPPUNIT_ASSERT_MESSAGE("test hasMemory == false", this->smxml->hasMemory() == false);
  delete smxml;
  cout << "SiconosMemoryXMLTest::testHasMemory - test2 .........................OK" << endl;
}

