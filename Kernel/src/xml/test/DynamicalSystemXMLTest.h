#ifndef __DynamicalSystemXMLTest__
#define __DynamicalSystemXMLTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "DynamicalSystemXML.h"
#include<iostream>
#include<vector>


class DynamicalSystemXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  DynamicalSystemXML ds;
  SiconosMatrix matrixRef;
  /*SiconosVector*/
  SimpleVector vectorRef;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(DynamicalSystemXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer
  CPPUNIT_TEST(testGetNumber);
  CPPUNIT_TEST(testGetType);
  CPPUNIT_TEST(testGetId);
  CPPUNIT_TEST(testGetN);
  CPPUNIT_TEST(testGetStepsInMemory);
  CPPUNIT_TEST(testGetX);
  CPPUNIT_TEST(testGetPluginName);

  //CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testGetNumber();
  void testGetType();
  void testGetId();
  void testGetN();
  void testGetStepsInMemory();
  void testGetX();
  void testGetPluginName();

public:
  void setUp();
  void tearDown();

};

#endif
