//$id$

#include "SiconosMemoryXMLTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryXMLTest);


void SiconosMemoryXMLTest::setUp()
{

  try
  {
    this->doc = xmlParseFile("SiconosMemory.xml");
    this->root = xmlDocGetRootElement(doc);
    //    this->smxml = SiconosMemoryXML(root);
  }
  catch (SiconosException e)
  {
    cout << "Error in SiconosMemoryXMLTest : " << e.report() << endl;
    exit(0);
  }

}

void SiconosMemoryXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}


//______________________________________________________________________________

void SiconosMemoryXMLTest::testHasMemory()
{
  this->smxml = new SiconosMemoryXML(SiconosDOMTreeTools::findNodeChild(this->root, "test1"));
  CPPUNIT_ASSERT_MESSAGE("test hasMemory == true", this->smxml->hasMemory() == true);
  delete smxml;
  cout << "SiconosMemoryXMLTest::testHasMemory - test1 .........................OK" << endl;


  this->smxml = new SiconosMemoryXML(SiconosDOMTreeTools::findNodeChild(this->root, "test2"));
  CPPUNIT_ASSERT_MESSAGE("test hasMemory == false", this->smxml->hasMemory() == false);
  delete smxml;
  cout << "SiconosMemoryXMLTest::testHasMemory - test2 .........................OK" << endl;
}

