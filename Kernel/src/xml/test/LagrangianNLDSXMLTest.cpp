#include "LagrangianNLDSXMLTest.h"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianNLDSXMLTest);


void LagrangianNLDSXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("DynamicalSystem.xml");
    this->root = xmlDocGetRootElement(doc);
    this->ds = LagrangianNLDSXML(root, false);

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
    cout << "Error in LagrangianNLDSXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void LagrangianNLDSXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void LagrangianNLDSXMLTest::testGetNdof()
{
  CPPUNIT_ASSERT_MESSAGE("testGetNdof : ", ds.getNdof() == 3);
  cout << "LagrangianNLDSXMLTest >>> testGetNdof ............................... OK\n ";
}

void LagrangianNLDSXMLTest::testGetQVelocity()
{
  CPPUNIT_ASSERT_MESSAGE("testGetQ : ", ds.getQ() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetQ0 : ", ds.getQ0() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetVelocity : ", ds.getVelocity() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetVelocity0 : ", ds.getVelocity0() == vectorRef);
  cout << "LagrangianNLDSXMLTest >>> testGetQVelocity .......................... OK\n ";
}


void LagrangianNLDSXMLTest::testGetMemory()
{
  cout << "LagrangianNLDSXMLTest >>> testGetMemory ***** MUST BE REIMPLEMENTED WITH THE NEW OBJECT SICONOSMEMORY !\n ";
}

void LagrangianNLDSXMLTest::testGetMass()
{
  SiconosMatrix matrixRef("matrix.dat", true);
  CPPUNIT_ASSERT_MESSAGE("test isMPlugin : ", ds.getMMatrix() == matrixRef);
  cout << "LagrangianNLDSXMLTest >>> testGetMass ............................... OK\n ";
}

void LagrangianNLDSXMLTest::testIsPlugin()
{
  CPPUNIT_ASSERT_MESSAGE("test isMPlugin : ", ds.isMPlugin() == false);
  CPPUNIT_ASSERT_MESSAGE("test isQNLInertiaPlugin : ", ds.isQNLInertiaPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isFintPlugin : ", ds.isFintPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isFextPlugin : ", ds.isFextPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianQFintPlugin : ", ds.isJacobianQFintPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianVelocityFintPlugin : ", ds.isJacobianVelocityFintPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianQQNLInertiaPlugin : ", ds.isJacobianQQNLInertiaPlugin());
  CPPUNIT_ASSERT_MESSAGE("test isJacobianVelocityQNLInertiaPlugin : ", ds.isJacobianVelocityQNLInertiaPlugin());

  cout << "LagrangianNLDSXMLTest >>> testIsPlugin .............................. OK\n ";
}

void LagrangianNLDSXMLTest::testGetPluginName()
{
  CPPUNIT_ASSERT_MESSAGE("test getQNLInertiaPlugin : ", ds.getQNLInertiaPlugin() == "BasicPlugin:computeQNLInertia");
  CPPUNIT_ASSERT_MESSAGE("test getFintPlugin : ", ds.getFintPlugin() == "BasicPlugin:computeFInt");
  CPPUNIT_ASSERT_MESSAGE("test getFextPlugin : ", ds.getFextPlugin() == "BasicPlugin:computeFExt");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianQFintPlugin : ", ds.getJacobianQFintPlugin() == "BasicPlugin:computeJacobianQFInt");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianVelocityFintPlugin : ", ds.getJacobianVelocityFintPlugin() == "BasicPlugin:computeJacobianVelocityFInt");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianQQNLInertiaPlugin : ", ds.getJacobianQQNLInertiaPlugin() == "BasicPlugin:computeJacobianQQNLInertia");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianVelocityQNLInertiaPlugin : ", ds.getJacobianVelocityQNLInertiaPlugin() == "BasicPlugin:computeJacobianVelocityQNLInertia");

  cout << "LagrangianNLDSXMLTest >>> testGetPluginName ......................... OK\n ";
}

