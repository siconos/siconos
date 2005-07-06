#include "DynamicalSystemXMLTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(DynamicalSystemXMLTest);




void DynamicalSystemXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("DynamicalSystem.xml");
    this->root = xmlDocGetRootElement(doc);
    this->ds = DSXML(root, false);

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
    cout << "Error in DynamicalSystemXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void DynamicalSystemXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void DynamicalSystemXMLTest::testGetNumber()
{
  CPPUNIT_ASSERT_MESSAGE("testGetNumber : ", ds.getNumber() == 120);
  cout << "DynamicalSystemXMLTest >>> testGetNumber ............................ OK\n ";
}

void DynamicalSystemXMLTest::testGetType()
{
  CPPUNIT_ASSERT_MESSAGE("testGetType : ", ds.getType() == "LagrangianDS");
  cout << "DynamicalSystemXMLTest >>> testGetType .............................. OK\n ";
}

void DynamicalSystemXMLTest::testGetId()
{
  CPPUNIT_ASSERT_MESSAGE("testGetId : ", ds.getId() == "Ball");
  cout << "DynamicalSystemXMLTest >>> testGetId ................................ OK\n ";
}

void DynamicalSystemXMLTest::testGetN()
{
  CPPUNIT_ASSERT_MESSAGE("testGetN : ", ds.getN() == 6);
  cout << "DynamicalSystemXMLTest >>> testGetN ................................. OK\n ";
}

void DynamicalSystemXMLTest::testGetStepsInMemory()
{
  CPPUNIT_ASSERT_MESSAGE("testGetStepsInMemory : ", ds.getStepsInMemory() == 3);
  cout << "DynamicalSystemXMLTest >>> testGetStepsInMemory ..................... OK\n ";
}

void DynamicalSystemXMLTest::testGetX()
{

  CPPUNIT_ASSERT_MESSAGE("testGetX : ", ds.getX() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetX0 : ", ds.getX0() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetXDot : ", ds.getXDot() == vectorRef);

  //  vector<SiconosVector*> XMem;
  //  XMem = ds.getXMemory();
  //
  //  CPPUNIT_ASSERT_MESSAGE("testGetXDot : ", *XMem[1]==vectorRef);

  cout << "DynamicalSystemXMLTest >>> testGetX ................................. OK\n ";
}

void DynamicalSystemXMLTest::testGetPluginName()
{
  CPPUNIT_ASSERT_MESSAGE("test getVectorFieldPlugin : ", ds.getVectorFieldPlugin() == "BasicPlugin:vectorField");
  CPPUNIT_ASSERT_MESSAGE("test getComputeJacobianXPlugin : ", ds.getComputeJacobianXPlugin() == "BasicPlugin:computeJacobianX");
  cout << "DynamicalSystemXMLTest >>> testGetPluginName ........................ OK\n ";
}



