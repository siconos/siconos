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


