#ifndef __TimeDiscretisationXMLTest__
#define __TimeDiscretisationXMLTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "TimeDiscretisationXML.h"


class TimeDiscretisationXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  TimeDiscretisationXML td;
  /*SiconosVector*/
  SimpleVector vectorRef;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(TimeDiscretisationXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer
  CPPUNIT_TEST(testIsConstant);
  CPPUNIT_TEST(testH);
  CPPUNIT_TEST(testgetN);

  CPPUNIT_TEST(testGetTk);
  CPPUNIT_TEST(testGetHminHmax);
  //CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testIsConstant();
  void testH();
  void testgetN();
  void testGetTk();
  void testGetHminHmax();


public:
  void setUp();
  void tearDown();

};

#endif
