#ifndef __RELAYNSLAWXMLTEST__
#define __RELAYNSLAWXMLTEST__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "RelayNSLXML.h"
#include "XMLException.h"




class RelayNSLXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  RelayNSLXML RNonSmoothLaw;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(RelayNSLXMLTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer

  CPPUNIT_TEST(testGetC);
  CPPUNIT_TEST(testGetD);
  CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);


  // on termine
  CPPUNIT_TEST_SUITE_END();



  // declaration de fonctions de test
  void testGetC();
  void testGetD();
  void testIfTagIsNotPresent();


public:
  void setUp();
  void tearDown();



};

#endif
