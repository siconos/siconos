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
#include "TimeDiscretisationXMLTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(TimeDiscretisationXMLTest);




void TimeDiscretisationXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("TimeDiscretisation.xml");
    this->root = xmlDocGetRootElement(doc);
    this->td = TimeDiscretisationXML(root);
  }
  catch (SiconosException e)
  {
    cout << "Error in TimeDiscretisationXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void TimeDiscretisationXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void TimeDiscretisationXMLTest::testIsConstant()
{
  CPPUNIT_ASSERT_MESSAGE("test isConstant", td.isConstant());
  cout << "TimeDiscretisationXMLTest >>> testIsConstant ................................. OK\n ";
}

void TimeDiscretisationXMLTest::testH()
{
  CPPUNIT_ASSERT_MESSAGE("test H value", td.getH() == 0.1);

  cout << "TimeDiscretisationXMLTest >>> testH ................................. OK\n ";
}

void TimeDiscretisationXMLTest::testgetN()
{
  CPPUNIT_ASSERT_MESSAGE("test H value", td.getN() == 1);

  cout << "TimeDiscretisationXMLTest >>> testGetN .............................. OK\n ";
}

void TimeDiscretisationXMLTest::testGetTk()
{
  vector<double> v(2);
  v.at(0) = 5;
  v.at(1) = -5;
  vectorRef = /*SiconosVector*/SimpleVector(v);

  //cout<<"td.getTk() : "<< *(td.getTk()) <<endl<<" - vectorRef : "<< vectorRef <<endl;

  CPPUNIT_ASSERT_MESSAGE("test tk value", td.getTk() == vectorRef);

  cout << "TimeDiscretisationXMLTest >>> testGetTk ............................. OK\n ";
}

void TimeDiscretisationXMLTest::testGetHminHmax()
{
  CPPUNIT_ASSERT_MESSAGE("test HMin value", td.getHMin() == 0);
  CPPUNIT_ASSERT_MESSAGE("test HMin value", td.getHMax() == 5);

  cout << "TimeDiscretisationXMLTest >>> testGetHminHmax ....................... OK\n ";
}


