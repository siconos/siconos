//$Id $
#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "CompositeVector.h"
#include "SiconosMatrix.h"

#include <math.h>
//#include "SiconosMatrix.h"

using namespace std;

class CompositeVectorTest : public CppUnit::TestFixture
{


private:

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(CompositeVectorTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testBuildCompositeVector);
  CPPUNIT_TEST(testBuildCompositeVector1);
  CPPUNIT_TEST(testOperatorAccess);
  CPPUNIT_TEST(testAdd);
  CPPUNIT_TEST(testSize);
  CPPUNIT_TEST(testWrite);
  CPPUNIT_TEST(testRead);

  CPPUNIT_TEST(testOperatorPlusEqualGEN);
  CPPUNIT_TEST(testOperatorSubEqualGEN);
  CPPUNIT_TEST(testOperatorEqualGEN);
  CPPUNIT_TEST(testOperatorComp);
  CPPUNIT_TEST(testOperatorCompDiff);
  CPPUNIT_TEST(testOperatorPlusEqualSPC);
  CPPUNIT_TEST(testOperatorSubEqualSPC);
  CPPUNIT_TEST(testOperatorMultEqualSPC);
  CPPUNIT_TEST(testOperatorDivEqualSPC);
  CPPUNIT_TEST(testOperatorEqualSPC);
  CPPUNIT_TEST(testExternalOperatorPlusGEN);
  CPPUNIT_TEST(testExternalOperatorSubGEN);
  CPPUNIT_TEST(testExternalOperatorMultDoubleSPC);
  CPPUNIT_TEST(testExternalOperatorDivDoubleSPC);
  CPPUNIT_TEST(testExternalOperatorPlusSPC);
  CPPUNIT_TEST(testExternalOperatorSubSPC);
  CPPUNIT_TEST(testExternalOperatorMultMat);
  CPPUNIT_TEST(testExternalOperatorMultTransMat);



  // les tests qui doivent echouer
  //CPPUNIT_TEST_FAIL(testFail);

  // on termine
  CPPUNIT_TEST_SUITE_END();



  // declaration de fonctions de test
  void testBuildCompositeVector();
  void testBuildCompositeVector1();
  void testOperatorAccess();
  void testAdd();
  void testSize();
  void testWrite();
  void testRead();

  void testOperatorPlusEqualGEN();
  void testOperatorSubEqualGEN();
  void testOperatorEqualGEN();
  void testOperatorComp();
  void testOperatorCompDiff();
  void testOperatorPlusEqualSPC();
  void testOperatorSubEqualSPC();
  void testOperatorMultEqualSPC();
  void testOperatorDivEqualSPC();
  void testOperatorEqualSPC();
  void testExternalOperatorPlusGEN();
  void testExternalOperatorSubGEN();
  void testExternalOperatorMultDoubleSPC();
  void testExternalOperatorDivDoubleSPC();
  void testExternalOperatorPlusSPC();
  void testExternalOperatorSubSPC();
  void testExternalOperatorMultMat();
  void testExternalOperatorMultTransMat();

  // declararion des variables de tests
  SimpleVector q;
  SimpleVector dotq;

public:
  void setUp();
  void tearDown();

};

#endif

//$Log $