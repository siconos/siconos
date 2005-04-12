//$id$

#ifndef PLATFORMTEST_H
#define PLATFORMTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "Model.h"
#include "LCP.h"
#include "RuntimeException.h"
#include "XMLException.h"

using namespace std;

class PlatformTest : public CppUnit::TestFixture
{
private:
  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(PlatformTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testOptionalAttributes1);
  CPPUNIT_TEST(testOptionalAttributes2);
  CPPUNIT_TEST(testOptionalAttributes3);
  CPPUNIT_TEST(testOptionalAttributes4);
  CPPUNIT_TEST(testOptionalAttributes5);
  CPPUNIT_TEST(testOptionalAttributes6);
  CPPUNIT_TEST(testOptionalAttributes7);
  CPPUNIT_TEST(testOptionalAttributes8);
  CPPUNIT_TEST(testOptionalAttributes9);
  CPPUNIT_TEST(testOptionalAttributes10);
  CPPUNIT_TEST(testOptionalAttributes11);
  CPPUNIT_TEST(testOptionalAttributes12);
  CPPUNIT_TEST(testOptionalAttributes13);
  CPPUNIT_TEST(testOptionalAttributes14);
  CPPUNIT_TEST(testOptionalAttributes15);
  CPPUNIT_TEST(testOptionalAttributes16);
  CPPUNIT_TEST(testManualCreation);
  CPPUNIT_TEST(testManualCreation2);
  CPPUNIT_TEST(testManualCreation3);
  CPPUNIT_TEST(testManualCreation4);
  CPPUNIT_TEST(testCheckXMLPlatform);
  CPPUNIT_TEST(testMixteCreation);
  CPPUNIT_TEST(testXMLSchemaAttributeGood);
  CPPUNIT_TEST(testXMLSchemaAttributeGood2);
  CPPUNIT_TEST(testXMLSchemaAttributeBad1);
  CPPUNIT_TEST_EXCEPTION(testXMLSchemaAttributeBad2, XMLException);
  CPPUNIT_TEST_EXCEPTION(testPlatformException, RuntimeException);
  CPPUNIT_TEST(testOSNSP1);
  CPPUNIT_TEST(testOSNSP2);
  CPPUNIT_TEST(testOSNSP3);
  CPPUNIT_TEST(testOSNSP4);
  CPPUNIT_TEST(testOSNSP5);
  CPPUNIT_TEST(testLmgc90);
  //CPPUNIT_TEST( testMainSiconos );


  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testOptionalAttributes1();
  void testOptionalAttributes2();
  void testOptionalAttributes3();
  void testOptionalAttributes4();
  void testOptionalAttributes5();
  void testOptionalAttributes6();
  void testOptionalAttributes7();
  void testOptionalAttributes8();
  void testOptionalAttributes9();
  void testOptionalAttributes10();
  void testOptionalAttributes11();
  void testOptionalAttributes12();
  void testOptionalAttributes13();
  void testOptionalAttributes14();
  void testOptionalAttributes15();
  void testOptionalAttributes16();
  void testManualCreation();
  void testManualCreation2();
  void testManualCreation3();
  void testManualCreation4();
  void testCheckXMLPlatform();
  void testMixteCreation();
  void testXMLSchemaAttributeGood();
  void testXMLSchemaAttributeGood2();
  void testXMLSchemaAttributeBad1();
  void testXMLSchemaAttributeBad2();
  void testPlatformException();
  void testOSNSP1();
  void testOSNSP2();
  void testOSNSP3();
  void testOSNSP4();
  void testOSNSP5();
  void testLmgc90();

  void testMainSiconos();
  //void testFail();

  // declararion des variables de tests


public:
  void setUp();
  void tearDown();
};

#endif // PLATFORMTEST_H

