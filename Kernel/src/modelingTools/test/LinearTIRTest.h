#ifndef __LinearTIRTest__
#define __LinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LinearTIR.h"

class LinearTIRTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LinearTIRTest);

  // tests to be done ...

  //CPPUNIT_TEST(testBuildLinearTIR);
  CPPUNIT_TEST(testBuildLinearTIR1);
  CPPUNIT_TEST(testBuildLinearTIR2);
  CPPUNIT_TEST(testSetC);
  CPPUNIT_TEST(testSetCPtr);
  CPPUNIT_TEST(testSetD);
  CPPUNIT_TEST(testSetDPtr);
  CPPUNIT_TEST(testSetF);
  CPPUNIT_TEST(testSetFPtr);
  CPPUNIT_TEST(testSetE);
  CPPUNIT_TEST(testSetEPtr);
  CPPUNIT_TEST(testSetB);
  CPPUNIT_TEST(testSetBPtr);
  CPPUNIT_TEST(testSetA);
  CPPUNIT_TEST(testSetAPtr);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLinearTIR();
  void testBuildLinearTIR1();
  void testBuildLinearTIR2();
  void testSetC();
  void testSetCPtr();
  void testSetD();
  void testSetDPtr();
  void testSetF();
  void testSetFPtr();
  void testSetE();
  void testSetEPtr();
  void testSetB();
  void testSetBPtr();
  void testSetA();
  void testSetAPtr();

  // Members

  SiconosMatrix *C, *B, *F, *D;
  SimpleVector *a, *e;

public:
  void setUp();
  void tearDown();

};

#endif




