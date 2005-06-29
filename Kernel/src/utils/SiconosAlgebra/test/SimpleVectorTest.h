#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "CompositeVector.h"
#include <math.h>
#include <vector>

class SimpleVectorTest : public CppUnit::TestFixture
{


private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SimpleVectorTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildSimpleVector);
  CPPUNIT_TEST(testBuildSimpleVector1);
  CPPUNIT_TEST(testBuildSimpleVector2);
  CPPUNIT_TEST(testBuildSimpleVector3);
  CPPUNIT_TEST(testBuildSimpleVector4);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testSetValue);
  CPPUNIT_TEST(testGetValue);
  CPPUNIT_TEST(testSetValues);
  CPPUNIT_TEST(testGetValues);
  CPPUNIT_TEST(testSize);
  CPPUNIT_TEST(testReadWrite);
  CPPUNIT_TEST(testOperatorPlusEqual);
  CPPUNIT_TEST(testOperatorEqual);
  CPPUNIT_TEST(testOperatorComp);
  CPPUNIT_TEST(testOperatorMultDivEqual);
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testSubtraction);
  CPPUNIT_TEST(testExternalOperatorPlusMoins);
  CPPUNIT_TEST(testExternalOperatorMultDiv);
  CPPUNIT_TEST(testExternalOperatorMultMat);
  CPPUNIT_TEST(testExternalOperatorMatTransMult);
  CPPUNIT_TEST(testOperatorAccess);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildSimpleVector();
  void testBuildSimpleVector1();
  void testBuildSimpleVector2();
  void testBuildSimpleVector3();
  void testBuildSimpleVector4();
  void testZero();
  void testSetValue();
  void testGetValue();
  void testSetValues();
  void testGetValues();
  void testSize();
  void testReadWrite();
  void testOperatorPlusEqual();
  void testOperatorEqual();
  void testOperatorComp();
  void testOperatorMultDivEqual();
  void testAddition();
  void testSubtraction();
  void testExternalOperatorPlusMoins();
  void testExternalOperatorMultDiv();
  void testExternalOperatorMultMat();
  void testExternalOperatorMatTransMult();
  void testOperatorAccess();

  // Members

  std::vector<double> vq;
  std::vector<double> vdotq;
  SimpleVector * u, *u2, *u3, *u4;
  SiconosVector *nsv;

public:
  void setUp();
  void tearDown();

};

#endif




