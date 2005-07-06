#ifndef LINEARSYSTEMDSXML_H
#define LINEARSYSTEMDSXML_H

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include "LinearDSXML.h"

using namespace std;

class LinearsystemDSXMLTest: public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root, *child;
  SiconosMatrix matrixRef;
  /*SiconosVector*/
  SimpleVector vectorRef;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(LinearsystemDSXMLTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testGetAB);
  CPPUNIT_TEST(testGetUF);

  // exceptions
  //CPPUNIT_TEST_EXCEPTION(testGetStringAttributeValueException, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();


  // declaration de fonctions de test
  void testGetAB();
  void testGetUF();

public:
  void setUp();
  void tearDown();

};

#endif // LINEARSYSTEMDSXML_H
