#include "LinearDSXMLTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(LinearsystemDSXMLTest);




void LinearsystemDSXMLTest::setUp()
{
  doc = xmlParseFile("LinearDS.xml");
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
  LinearDSXML lsdsxml(root, false);
  SiconosMatrix A, E;
  A = lsdsxml.getA();
  E = lsdsxml.getE();

  CPPUNIT_ASSERT_MESSAGE("testGetAB : A == matrixRef", A == matrixRef);
  CPPUNIT_ASSERT_MESSAGE("testGetAB : E == matrixRef", E == matrixRef);
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
  LinearDSXML lsdsxml(root, false);
  /*SiconosVector*/
  SimpleVector U, b;
  U = lsdsxml.getUVector();
  b = lsdsxml.getBVector();
  CPPUNIT_ASSERT_MESSAGE("testGetUF : U == vectorRef", U == vectorRef);
  CPPUNIT_ASSERT_MESSAGE("testGetUF : b == vectorRef", b == vectorRef);
  cout << " LinearsystemDSXMLTest >>> testGetUF ................................ OK\n ";
}
