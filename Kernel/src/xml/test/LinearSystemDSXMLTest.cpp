#include "LinearSystemDSXMLTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(LinearsystemDSXMLTest);




void LinearsystemDSXMLTest::setUp()
{
  doc = xmlParseFile("LinearSystemDS.xml");
  root = xmlDocGetRootElement(doc);
  child = NULL;
  matrixRef = SiconosMatrix("matrix.dat", true);
  vectorRef = /*SiconosVector*/SimpleVector("vector.dat", true);
}

void LinearsystemDSXMLTest::tearDown()
{
  //xmlFreeDoc(doc);
  xmlCleanupParser();
}

void LinearsystemDSXMLTest::testGetAB()
{
  //  try{
  LinearSystemDSXML lsdsxml(root, false);
  SiconosMatrix A, B;
  A = lsdsxml.getA();
  B = lsdsxml.getB();

  CPPUNIT_ASSERT_MESSAGE("testGetAB : A == matrixRef", A == matrixRef);
  CPPUNIT_ASSERT_MESSAGE("testGetAB : B == matrixRef", B == matrixRef);
  //  }
  //  catch(SiconosException e)
  //  {
  //    cout << "LinearsystemDSXMLTest error : "<<e.report() <<endl;
  //    exit(0);
  //  }

  cout << " LinearsystemDSXMLTest >>> testGetAB ................................ OK\n ";
}

void LinearsystemDSXMLTest::testGetUF()
{
  LinearSystemDSXML lsdsxml(root, false);
  /*SiconosVector*/
  SimpleVector U, F;
  U = lsdsxml.getUVector();
  F = lsdsxml.getFVector();
  CPPUNIT_ASSERT_MESSAGE("testGetUF : U == vectorRef", U == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetUF : F == vectorRef", F == vectorRef);
  cout << " LinearsystemDSXMLTest >>> testGetUF ................................ OK\n ";
}
