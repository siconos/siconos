//$Id: SimpleVectorTest.h,v 1.11 2004/09/14 13:24:54 charlety Exp $
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

//$Log: SimpleVectorTest.h,v $
//Revision 1.11  2004/09/14 13:24:54  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.10  2004/09/09 14:32:50  charlety
//
//_ New tests for operators of multiplication between vectors and matrices.
//
//Revision 1.9  2004/08/24 11:29:20  charlety
//
//_ methods replacing generic operators for mixed operations done.
//
//Revision 1.8  2004/08/23 09:29:02  charlety
//
//_ tests for compositeVector in progress
//
//Revision 1.7  2004/08/20 15:00:35  charlety
//
//_ Tests for operators of SimpleVector
//
//Revision 1.6  2004/08/19 15:21:28  charlety
//
//_ SimpleVector and CompositeVector in progress.
//_ for the operators, we prefer now using directly functions of Blas1++ instead
//  of these of Blas++.h
//
//Revision 1.5  2004/08/16 09:40:01  charlety
//
//_ new tests for the simpleVector
//
//Revision 1.4  2004/08/13 10:36:11  charlety
//
//_ tests for simpleVector in progress
//
//Revision 1.3  2004/08/11 14:16:08  charlety
//
//_ NewSiconosVector in progress...(NewSiconosVector is an abstract class and
//  SimpleVector inherits of NewSiconosVector).
//
//Revision 1.2  2004/07/30 14:21:55  charlety
//
//_ new functions and tests for the new SiconosVector
//