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
#include "RelayNSLXMLTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(RelayNSLXMLTest);




void RelayNSLXMLTest::setUp()
{
  this->doc = xmlParseFile("RelayNSL.xml");
  this->root = xmlDocGetRootElement(doc);
  this->RNonSmoothLaw = RelayNSLXML(root);
}

void RelayNSLXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void RelayNSLXMLTest::testGetC()
{
  double c = RNonSmoothLaw.getC();

  CPPUNIT_ASSERT_MESSAGE("testGetC : c", c == 0.054);
  cout << "RelayNSLXMLTest >>> testGetC ...................................... OK\n ";
}

void RelayNSLXMLTest::testGetD()
{
  double d = RNonSmoothLaw.getD();

  CPPUNIT_ASSERT_MESSAGE("testGetD : d", d == 0.064);
  cout << "RelayNSLXMLTest >>> testGetD ...................................... OK\n ";
}


void RelayNSLXMLTest::testIfTagIsNotPresent()
{
  xmlFreeDoc(doc);
  xmlCleanupParser();
  this->doc = xmlParseFile("RelayNSLCTagIsNotPresent.xml");
  this->root = xmlDocGetRootElement(doc);
  this->RNonSmoothLaw = RelayNSLXML(root);


}


