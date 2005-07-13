#ifndef __NonSmoothDynamicalSystemTest__
#define __NonSmoothDynamicalSystemTest__

#include <cppunit/extensions/HelperMacros.h>
#include "NonSmoothDynamicalSystem.h"
#include "LagrangianLinearTIDS.h"

class NonSmoothDynamicalSystemTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(NonSmoothDynamicalSystemTest);

  // tests to be done ...

  //CPPUNIT_TEST(testBuildNonSmoothDynamicalSystem);
  CPPUNIT_TEST(testBuildNonSmoothDynamicalSystem1);
  CPPUNIT_TEST(testBuildNonSmoothDynamicalSystem2);
  CPPUNIT_TEST(testaddDynamicalSystem);
  CPPUNIT_TEST(testaddInteraction);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildNonSmoothDynamicalSystem1();
  void testBuildNonSmoothDynamicalSystem2();
  void testaddDynamicalSystem();
  void testaddInteraction();

  // Members

  xmlNode * node;
  NSDSXML* tmpxml;
public:
  void setUp();
  void tearDown();

};

#endif




