#ifndef INTERACTIONXMLTEST_H
#define INTERACTIONXMLTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include "InteractionXML.h"

using namespace std;

class InteractionXMLTest: public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root, *child;
  InteractionXML interaction;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(InteractionXMLTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testGetNumber);
  CPPUNIT_TEST(testGetStatus);
  CPPUNIT_TEST(testGetId);
  CPPUNIT_TEST(testHasYLambda);
  CPPUNIT_TEST(testGetDSConcerned);

  // exceptions
  //CPPUNIT_TEST_EXCEPTION(testGetStringAttributeValueException, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();



  // declaration de fonctions de test
  void testGetNumber();
  void testGetStatus();
  void testGetId();
  void testHasYLambda();
  void testGetDSConcerned();

public:
  void setUp();
  void tearDown();

};

#endif // INTERACTIONXMLTEST_H
