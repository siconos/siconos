#ifndef __OneStepIntegratorXMLTest__
#define __OneStepIntegratorXMLTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include "OneStepIntegratorXML.h"
#include "SiconosModelXML.h"

#include "AdamsXML.h"
#include "MoreauXML.h"
#include "LsodarXML.h"

using namespace std;

class OneStepIntegratorXMLTest: public CppUnit::TestFixture
{


private:
  SiconosModelXML* modelXML;
  StrategyXML* strategyXML;
  vector<OneStepIntegratorXML*> oneStepIs;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(OneStepIntegratorXMLTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testNbOneStepIntegrator);
  CPPUNIT_TEST(testGetR);
  CPPUNIT_TEST(testGetDSConcerned);
  CPPUNIT_TEST(testGetType);
  CPPUNIT_TEST(testAdamsXML);
  CPPUNIT_TEST(testMoreauXML);
  CPPUNIT_TEST(testLsodarXML);
  // exceptions
  //CPPUNIT_TEST_EXCEPTION(testGetStringAttributeValueException, XMLException);


  // on termine
  CPPUNIT_TEST_SUITE_END();


  // declaration de fonctions de testvoid OneStepIntegratorXMLTest::testGetR()
  void testNbOneStepIntegrator();
  void testGetR();
  void testGetDSConcerned();
  void testGetType();
  void testAdamsXML();
  void testMoreauXML();
  void testLsodarXML();


public:
  void setUp();
  void tearDown();

};

#endif // __OneStepIntegratorXMLTest__
