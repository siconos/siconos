#include "LagrangianDSXMLTest.h"


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

void LagrangianDSXMLTest::testGetPluginName()
{
  CPPUNIT_ASSERT_MESSAGE("test getQNLInertiaPlugin : ", ds.getQNLInertiaPlugin() == "BasicPlugin:computeQNLInertia");
  CPPUNIT_ASSERT_MESSAGE("test getFintPlugin : ", ds.getFintPlugin() == "BasicPlugin:computeFInt");
  CPPUNIT_ASSERT_MESSAGE("test getFextPlugin : ", ds.getFextPlugin() == "BasicPlugin:computeFExt");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianQFintPlugin : ", ds.getJacobianQFintPlugin() == "BasicPlugin:computeJacobianQFInt");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianVelocityFintPlugin : ", ds.getJacobianVelocityFintPlugin() == "BasicPlugin:computeJacobianVelocityFInt");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianQQNLInertiaPlugin : ", ds.getJacobianQQNLInertiaPlugin() == "BasicPlugin:computeJacobianQQNLInertia");
  CPPUNIT_ASSERT_MESSAGE("test getJacobianVelocityQNLInertiaPlugin : ", ds.getJacobianVelocityQNLInertiaPlugin() == "BasicPlugin:computeJacobianVelocityQNLInertia");

  cout << "LagrangianDSXMLTest >>> testGetPluginName ......................... OK\n ";
}

