#include "BoundaryConditionXMLTest.h"
#include<iostream>
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(BoundaryConditionXMLTest);




void BoundaryConditionXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("BoundaryConditionLinear.xml");
    this->root = xmlDocGetRootElement(doc);
    this->bcLinear = LinearBCXML(SiconosDOMTreeTools::findNodeChild(root));

    this->doc = xmlParseFile("BoundaryConditionPeriodic.xml");
    this->root = xmlDocGetRootElement(doc);
    this->bcPeriodic = PeriodicBCXML(SiconosDOMTreeTools::findNodeChild(root));

    this->doc = xmlParseFile("BoundaryConditionNonLinear.xml");
    this->root = xmlDocGetRootElement(doc);
    this->bcNLinear = NLinearBCXML(SiconosDOMTreeTools::findNodeChild(root));


    matrixRef = SiconosMatrix("matrix.dat", true);
    vectorRef = /*SiconosVector*/SimpleVector("vector.dat", true);
  }
  catch (SiconosException e)
  {
    cout << "Error in BoundaryConditionXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void BoundaryConditionXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void BoundaryConditionXMLTest::testGetType()
{

  CPPUNIT_ASSERT_MESSAGE("testGetType : linear", bcLinear.getType() == "Linear");
  CPPUNIT_ASSERT_MESSAGE("testGetType : non-linear", bcNLinear.getType() == "NLinear");
  CPPUNIT_ASSERT_MESSAGE("testGetType : periodic", bcPeriodic.getType() == "Periodic");

  cout << "BoundaryConditionXMLTest >>> testGetType ............................ OK\n ";
}

void BoundaryConditionXMLTest::testGetOmega()
{

  CPPUNIT_ASSERT_MESSAGE("testGetOmega : Omega", bcLinear.getOmega() == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetOmega : Omega0", bcLinear.getOmega0() == matrixRef);
  CPPUNIT_ASSERT_MESSAGE("testGetOmega : OmegaT", bcLinear.getOmegaT() == matrixRef);

  cout << "BoundaryConditionXMLTest >>> testGetOmega ........................... OK\n ";
}


