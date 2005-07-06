#ifndef __LINEARTIRELATIONXMLTEST__
#define __LINEARTIRELATIONXMLTEST__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "LinearTIRXML.h"
#include "XMLException.h"
#include <iostream>

class LinearTIRXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  LinearTIRXML LinearTIR;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(LinearTIRXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer

  CPPUNIT_TEST(testGetC);
  CPPUNIT_TEST(testGetD);
  CPPUNIT_TEST(testGeta);

  CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);
  CPPUNIT_TEST_EXCEPTION(testIfAttributeNotPresent, XMLException);


  // on termine
  CPPUNIT_TEST_SUITE_END();



  // declaration de fonctions de test
  void testGetC();
  void testGetD();
  void testGeta();
  void testIfTagIsNotPresent();
  void testIfAttributeNotPresent();


public:
  void setUp();
  void tearDown();



};

#endif
