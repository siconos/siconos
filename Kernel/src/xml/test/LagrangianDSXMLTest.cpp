/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "LagrangianDSXMLTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianDSXMLTest);


void LagrangianDSXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("DynamicalSystem.xml");
    this->root = xmlDocGetRootElement(doc);
    this->ds = LagrangianDSXML(root, false);

    vector<double> v(6);
    v.at(0) = 1.0;
    v.at(1) = 0.0;
    v.at(2) = 0.0;
    v.at(3) = 0.0;
    v.at(4) = 0.0;
    v.at(5) = 0.0;
    vectorRef = /*SiconosVector*/SimpleVector(v);

    matrixRef = SiconosMatrix("matrix.dat", true);
    //vectorRef = SiconosVector("vector.dat", true);
  }
  catch (SiconosException e)
  {
    cout << "Error in LagrangianDSXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void LagrangianDSXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void LagrangianDSXMLTest::testGetNdof()
{
  CPPUNIT_ASSERT_MESSAGE("testGetNdof : ", ds.getNdof() == 3);
  cout << "LagrangianDSXMLTest >>> testGetNdof ............................... OK\n ";
}

void LagrangianDSXMLTest::testGetQVelocity()
{
  CPPUNIT_ASSERT_MESSAGE("testGetQ : ", ds.getQ() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetQ0 : ", ds.getQ0() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetVelocity : ", ds.getVelocity() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetVelocity0 : ", ds.getVelocity0() == vectorRef);
  cout << "LagrangianDSXMLTest >>> testGetQVelocity .......................... OK\n ";
}


void LagrangianDSXMLTest::testGetMemory()
{
  cout << "LagrangianDSXMLTest >>> testGetMemory ***** MUST BE REIMPLEMENTED WITH THE NEW OBJECT SICONOSMEMORY !\n ";
}

void LagrangianDSXMLTest::testGetMass()
{
  SiconosMatrix matrixRef("matrix.dat", true);
  CPPUNIT_ASSERT_MESSAGE("test isMPlugin : ", ds.getMMatrix() == matrixRef);
  cout << "LagrangianDSXMLTest >>> testGetMass ............................... OK\n ";
}

void LagrangianDSXMLTest::testIsPlugin()
{
  CPPUNIT_ASSERT_MESSAGE("test isMPlugin : ", ds.isMPlugin() == false);
  CPPUNIT_ASSERT_MESSAGE("test isQNLInertiaPlugin : ", ds.isQNLInertiaPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isFintPlugin : ", ds.isFintPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isFextPlugin : ", ds.isFextPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianQFintPlugin : ", ds.isJacobianQFintPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianVelocityFintPlugin : ", ds.isJacobianVelocityFintPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianQQNLInertiaPlugin : ", ds.isJacobianQQNLInertiaPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianVelocityQNLInertiaPlugin : ", ds.isJacobianVelocityQNLInertiaPlugin());

  cout << "LagrangianDSXMLTest >>> testIsPlugin .............................. OK\n ";
}


