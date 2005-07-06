#ifndef __SiconosDOMTreeToolsTest__
#define __SiconosDOMTreeToolsTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"

using namespace std;

class SiconosDOMTreeToolsTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root, *child;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(SiconosDOMTreeToolsTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testGetStringContentValue);
  CPPUNIT_TEST(testGetIntegerContentValue);
  CPPUNIT_TEST(testGetDoubleContentValue);

  CPPUNIT_TEST(testGetStringAttributeValue);
  CPPUNIT_TEST(testGetDoubleAttributeValue);
  CPPUNIT_TEST(testGetIntegerAttributeValue);
  CPPUNIT_TEST(testGetBooleanAttributeValue);

  CPPUNIT_TEST(testGetSiconosVectorValue);
  CPPUNIT_TEST(testGetSiconosMatrixValue);
  CPPUNIT_TEST(testGetSiconosMatrixValue2);
  //  CPPUNIT_TEST(testGetVectorMemoryValue);

  CPPUNIT_TEST(testSetStringContentValue);
  CPPUNIT_TEST(testSetIntegerContentValue);
  CPPUNIT_TEST(testSetDoubleContentValue);
  CPPUNIT_TEST(testSetStringAttributeValue);
  CPPUNIT_TEST(testSetDoubleAttributeValue);
  CPPUNIT_TEST(testSetSiconosVectorValue);
  CPPUNIT_TEST(testSetSiconosMatrixValue);
  CPPUNIT_TEST(testSetIntegerAttributeValue);
  CPPUNIT_TEST(testSetBooleanAttributeValue);
  //  CPPUNIT_TEST(testSetVectorMemoryValue);

  CPPUNIT_TEST(testSetSiconosVectorValue);
  CPPUNIT_TEST(testGetSiconosVectorValue);
  //  CPPUNIT_TEST(testSetSiconosVectorValue);

  //  CPPUNIT_TEST(testString2Vector);

  CPPUNIT_TEST(testCreateMatrixNode);
  CPPUNIT_TEST(testCreateVectorNode);
  CPPUNIT_TEST(testCreateIntegerNode);
  CPPUNIT_TEST(testCreateDoubleNode);
  CPPUNIT_TEST(testCreateStringNode);
  CPPUNIT_TEST(testCreateBooleanNode);

  //  CPPUNIT_TEST(testGetNodeChildrenNumber);


  // exceptions
  CPPUNIT_TEST_EXCEPTION(testGetStringAttributeValueException, XMLException);
  CPPUNIT_TEST_EXCEPTION(testGetIntegerAttributeValueException, XMLException);
  CPPUNIT_TEST_EXCEPTION(testGetDoubleAttributeValueException, XMLException);
  CPPUNIT_TEST_EXCEPTION(testGetBooleanAttributeValueException, XMLException);

  CPPUNIT_TEST_EXCEPTION(testSetStringAttributeValueException, XMLException);
  CPPUNIT_TEST_EXCEPTION(testSetIntegerAttributeValueException, XMLException);
  CPPUNIT_TEST_EXCEPTION(testSetDoubleAttributeValueException, XMLException);
  CPPUNIT_TEST_EXCEPTION(testSetBooleanAttributeValueException, XMLException);
  CPPUNIT_TEST_EXCEPTION(testGetVectorValueBadSizeException, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();



  // declaration de fonctions de test
  void testGetStringContentValue();
  void testGetIntegerContentValue();
  void testGetDoubleContentValue();

  void testGetStringAttributeValue();
  void testGetDoubleAttributeValue();
  void testGetIntegerAttributeValue();
  void testGetBooleanAttributeValue();
  void testGetSiconosVectorValue();
  void testGetSiconosMatrixValue();
  void testGetSiconosMatrixValue2();
  //  void testGetVectorMemoryValue();

  void testSetStringContentValue();
  void testSetIntegerContentValue();
  void testSetDoubleContentValue();
  void testSetStringAttributeValue();
  void testSetDoubleAttributeValue();
  void testSetIntegerAttributeValue();
  void testSetSiconosVectorValue();
  void testSetSiconosMatrixValue();
  void testSetBooleanAttributeValue();
  //  void testSetVectorMemoryValue();

  //  void testString2Vector();

  void testCreateMatrixNode();
  void testCreateVectorNode();
  void testCreateIntegerNode();
  void testCreateDoubleNode();
  void testCreateStringNode();
  void testCreateBooleanNode();

  void testGetNodeChildrenNumber();

  void testGetStringAttributeValueException();
  void testGetIntegerAttributeValueException();
  void testGetDoubleAttributeValueException();
  void testGetBooleanAttributeValueException();

  void testSetStringAttributeValueException();
  void testSetIntegerAttributeValueException();
  void testSetDoubleAttributeValueException();
  void testSetBooleanAttributeValueException();
  void testGetVectorValueBadSizeException();

public:


  void setUp();
  void tearDown();

};

#endif
