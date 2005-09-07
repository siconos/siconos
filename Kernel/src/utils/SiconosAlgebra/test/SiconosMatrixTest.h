#ifndef __SiconosMatrixTest__
#define __SiconosMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
//#include "SiconosVector.h"
#include "SiconosMatrix.h"

using namespace std;

class SiconosMatrixTest : public CppUnit::TestFixture
{


private:

  // test suite
  CPPUNIT_TEST_SUITE(SiconosMatrixTest);

  CPPUNIT_TEST(testConstructor1);
  CPPUNIT_TEST(testConstructor2);
  CPPUNIT_TEST(copyConstructor);
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

  CPPUNIT_TEST_SUITE_END();


  void testConstructor1();
  void testConstructor2();
  void copyConstructor();
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

  SiconosMatrix A, B, C;
  SimpleVector SV;
  LaVectorDouble LVD;

public:
  void setUp();
  void tearDown();

};

#endif
