#include "RelayNSLXMLTest.h"
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


