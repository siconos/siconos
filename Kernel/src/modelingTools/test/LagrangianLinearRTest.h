#ifndef __LagrangianLinearRTest__
#define __LagrangianLinearRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LagrangianLinearR.h"

class LagrangianLinearRTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianLinearRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianLinearR1);
  CPPUNIT_TEST(testBuildLagrangianLinearR2);
  CPPUNIT_TEST(testBuildLagrangianLinearR3);
  CPPUNIT_TEST(testSetH);
  CPPUNIT_TEST(testSetHPtr);
  CPPUNIT_TEST(testSetB);
  CPPUNIT_TEST(testSetBPtr);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianLinearR1();
  void testBuildLagrangianLinearR2();
  void testBuildLagrangianLinearR3();
  void testSetH();
  void testSetHPtr();
  void testSetB();
  void testSetBPtr();
  void End();

  // Members

  SiconosMatrix *H;
  SimpleVector *b;

public:
  void setUp();
  void tearDown();

};

#endif




