#ifndef __SiconosMatrixTest__
#define __SiconosMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
//#include "SiconosVector.h"
#include "SiconosMatrix.h"

using namespace std;

class SiconosMatrixTest : public CppUnit::TestFixture
{


private:

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(SiconosMatrixTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(testEquality);
  CPPUNIT_TEST(testAffectation);
  CPPUNIT_TEST(testReadWriteAscii);
  CPPUNIT_TEST(testReadWriteBinary);
  CPPUNIT_TEST(testLinearSolve);
  CPPUNIT_TEST(testAddRow);
  CPPUNIT_TEST(testOperator);
  CPPUNIT_TEST(testBlockMatrixCopy1);
  CPPUNIT_TEST(testBlockMatrixCopy2);
  CPPUNIT_TEST(testBlockMatrixCopy3);
  CPPUNIT_TEST(testBlockMatrixCopy4);
  CPPUNIT_TEST_EXCEPTION(testSizeException, SiconosMatrixException);
  CPPUNIT_TEST_EXCEPTION(testConstructorException, SiconosMatrixException);
  CPPUNIT_TEST_EXCEPTION(testBlockMatrixCopyException1, SiconosMatrixException);
  CPPUNIT_TEST_EXCEPTION(testBlockMatrixCopyException2, SiconosMatrixException);


  // les tests qui doivent echouer
  //CPPUNIT_TEST_FAIL(testFail);

  // on termine
  CPPUNIT_TEST_SUITE_END();


  // declaration de fonctions de test
  void testConstructor1();
  void testConstructor2();
  void testEquality();
  void testAffectation();
  void testReadWriteAscii();
  void testReadWriteBinary();
  void testLinearSolve();
  void testAddRow();
  void testOperator();
  void testSizeException();
  void testConstructorException();

  void testBlockMatrixCopy1();
  void testBlockMatrixCopy2();
  void testBlockMatrixCopy3();
  void testBlockMatrixCopy4();

  void testBlockMatrixCopyException1();
  void testBlockMatrixCopyException2();

  // declararion des variables de tests
  SiconosMatrix A, B, C;
  /*SiconosVector SV;*/
  SimpleVector SV;
  LaVectorDouble LVD;

public:
  void setUp();
  void tearDown();

};

#endif
