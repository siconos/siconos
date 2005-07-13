#include "ModelTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ModelTest);


void ModelTest::setUp()
{
  t0 = 0.0;
  T = 10.0;
}

void ModelTest::tearDown()
{}

// xml constructor
void ModelTest::testBuildModel1()
{
  char * xmlFile = "ModelXml_test.xml" ;
  Model * M = new Model(xmlFile);
  // check that xml link is not NULL
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1a : ", M->getSiconosModelXMLPtr() != NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1b : ", M->getT0() == 0.1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1c : ", M->getFinalT() == 10.0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1d : ", M->getTitle() == "tryMxml", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1e : ", M->getAuthor() == "SiconosTeam", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1f : ", M->getDate() == "Today", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1g : ", M->getXmlSchema() == "/share/SICONOS/SiconosModelSchema-V1.2.xsd", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1h : ", M->getDescription() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1i : ", M->getNonSmoothDynamicalSystemPtr() != NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildModel1j : ", M->getStrategyPtr() != NULL, true);
  delete M;
  cout << " Constructor Model 1 (xml) ok" << endl;
}

// data constructor (1)
void ModelTest::testBuildModel2()
{
  Model * M = new Model(t0, T, "tryMxml", "SiconosTeam", "none", "Today", "/share/SICONOS/SiconosModelSchema-V1.2.xsd");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2a : ", M->getSiconosModelXMLPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2b : ", M->getT0() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2c : ", M->getFinalT() == 10.0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2d : ", M->getTitle() == "tryMxml", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2e : ", M->getAuthor() == "SiconosTeam", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2f : ", M->getDate() == "Today", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2g : ", M->getXmlSchema() == "/share/SICONOS/SiconosModelSchema-V1.2.xsd", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2h : ", M->getDescription() == "none", true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2i : ", M->getNonSmoothDynamicalSystemPtr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMode2j : ", M->getStrategyPtr() == NULL, true);
  delete M;
  cout << " Constructor Model 2 ok" << endl;
}

// set Strategy

void ModelTest::testSetStrategyPtr()
{
  Strategy * newStrat = new Strategy();

  Model * M = new Model(t0, T);
  M->setStrategyPtr(newStrat);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test SetStrategyPtr : ", M->getStrategyPtr() == newStrat, true);
  delete M;
  cout << " test setStrategyPtr ok" << endl;
}

// setNonSmoothDynamicalSystem
void ModelTest::testsetNonSmoothDynamicalSystemPtr()
{
  NonSmoothDynamicalSystem * newNsds = new NonSmoothDynamicalSystem();
  Model * M = new Model(t0, T);
  M->setNonSmoothDynamicalSystemPtr(newNsds);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test setNonSmoothDynamicalSystemPtr : ", M->getNonSmoothDynamicalSystemPtr() == newNsds, true);
  delete M;
  cout << " test setNonSmoothDynamicalSystemPtr ok " << endl;
}

// setSiconosModelXML
void ModelTest::testsetSiconosModelXMLPtr()
{
  SiconosModelXML * newMxml = new SiconosModelXML();
  Model * M = new Model(t0, T);
  M->setSiconosModelXMLPtr(newMxml);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" test setSiconosModelXMLPtr: ", M->getSiconosModelXMLPtr() == newMxml, true);
  delete M;
  cout << " test setSiconosModelXMLPtr ok " << endl;
}
