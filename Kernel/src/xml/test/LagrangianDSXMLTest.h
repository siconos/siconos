#ifndef __LagrangianDSXMLTest__
#define __LagrangianDSXMLTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "LagrangianDSXML.h"




class LagrangianDSXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  LagrangianDSXML ds;
  SiconosMatrix matrixRef;
  /*SiconosVector*/
  SimpleVector vectorRef;


  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(LagrangianDSXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer
  CPPUNIT_TEST(testGetNdof);
  CPPUNIT_TEST(testGetQVelocity);
  CPPUNIT_TEST(testGetMass);
  CPPUNIT_TEST(testIsPlugin);
  CPPUNIT_TEST(testGetPluginName);
  CPPUNIT_TEST(testGetMemory);

  //CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testGetNdof();
  void testGetQVelocity();
  void testGetMemory();
  void testGetMass();
  void testIsPlugin();
  void testGetPluginName();

public:
  void setUp();
  void tearDown();
};

#endif
