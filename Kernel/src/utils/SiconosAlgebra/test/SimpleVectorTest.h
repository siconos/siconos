#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include <math.h>

#include <time.h>
#include <stdio.h>


//#include "SiconosMatrix.h"

using namespace std;

class SimpleVectorTest : public CppUnit::TestFixture
{


private:

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(SimpleVectorTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testBuildSimpleVector);
  CPPUNIT_TEST(testBuildSimpleVector1);
  CPPUNIT_TEST(testBuildSimpleVector2);
  CPPUNIT_TEST(testBuildSimpleVector3);
  CPPUNIT_TEST(testOperatorEqual);
  CPPUNIT_TEST(testOperatorComp);
  CPPUNIT_TEST(testOperatorCompDiff);
  CPPUNIT_TEST(testOperatorAccessRef);
  CPPUNIT_TEST(testOperatorAccess);
  CPPUNIT_TEST(testOperatorPlusEqualGEN);
  CPPUNIT_TEST(testOperatorSubEqualGEN);
  CPPUNIT_TEST(testOperatorEqualGEN);
  CPPUNIT_TEST(testOperatorPlusEqualSPC);
  CPPUNIT_TEST(testOperatorSubEqualSPC);
  CPPUNIT_TEST(testOperatorMultEqualSPC);
  CPPUNIT_TEST(testOperatorDivEqualSPC);
  CPPUNIT_TEST(testExternalOperatorPlusGEN);
  CPPUNIT_TEST(testExternalOperatorSubGEN);
  CPPUNIT_TEST(testExternalOperatorMultDoubleSPC);
  CPPUNIT_TEST(testExternalOperatorDivDoubleSPC);
  CPPUNIT_TEST(testExternalOperatorPlusSPC);
  CPPUNIT_TEST(testExternalOperatorSubSPC);
  CPPUNIT_TEST(testExternalOperatorMultMat);
  CPPUNIT_TEST(testExternalOperatorMultTransMat);

  CPPUNIT_TEST(testSetValues);
  CPPUNIT_TEST(testSize);
  CPPUNIT_TEST(testWrite);
  CPPUNIT_TEST(testRead);

  CPPUNIT_TEST(testDisplay);
  CPPUNIT_TEST(testMiscelleanous);
  //  CPPUNIT_TEST(testLapackVSblasMultDouble);


  CPPUNIT_TEST_EXCEPTION(testOperatorAccessRefException, SiconosVectorException);
  CPPUNIT_TEST_EXCEPTION(testOperatorAccessRefException1, SiconosVectorException);
  CPPUNIT_TEST_EXCEPTION(testOperatorAccessRefException2, SiconosVectorException);
  CPPUNIT_TEST_EXCEPTION(testOperatorAccessException, SiconosVectorException);
  CPPUNIT_TEST_EXCEPTION(testOperatorAccessException1, SiconosVectorException);
  CPPUNIT_TEST_EXCEPTION(testOperatorAccessException2, SiconosVectorException);
  //  CPPUNIT_TEST_EXCEPTION(testAddException, SiconosVectorException);// CPPUNIT_TEST_EXCEPTION(testFriendOperatorPlusException, SiconosVectorException);

  // les tests qui doivent echouer
  //CPPUNIT_TEST_FAIL(testFail);

  // on termine
  CPPUNIT_TEST_SUITE_END();



  // declaration de fonctions de test
  void testBuildSimpleVector();
  void testBuildSimpleVector1();
  void testBuildSimpleVector2();
  void testBuildSimpleVector3();
  void testOperatorAccessRef();
  void testOperatorAccess();
  void testSetValues();
  void testSize();
  void testWrite();
  void testRead();

  void testOperatorEqual();
  void testOperatorComp();
  void testOperatorCompDiff();
  void testOperatorPlusEqualGEN();
  void testOperatorSubEqualGEN();
  void testOperatorEqualGEN();
  void testOperatorPlusEqualSPC();
  void testOperatorSubEqualSPC();
  void testOperatorMultEqualSPC();
  void testOperatorDivEqualSPC();
  void testExternalOperatorPlusGEN();
  void testExternalOperatorSubGEN();
  void testExternalOperatorMultDoubleSPC();
  void testExternalOperatorDivDoubleSPC();
  void testExternalOperatorPlusSPC();
  void testExternalOperatorSubSPC();
  void testExternalOperatorMultMat();
  void testExternalOperatorMultTransMat();

  void testOperatorAccessRefException();
  void testOperatorAccessRefException1();
  void testOperatorAccessRefException2();
  void testOperatorAccessException();
  void testOperatorAccessException1();
  void testOperatorAccessException2();
  //  void testAddException();

  void testDisplay();
  void testMiscelleanous();
  void testLapackVSblasMultDouble();


  // declararion des variables de tests

public:
  void setUp();
  void tearDown();

};

#endif

