#ifndef __ModelTest__
#define __ModelTest__

#include <cppunit/extensions/HelperMacros.h>
#include "Model.h"
#include "TimeStepping.h"
#include "EventDriven.h"

class ModelTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(ModelTest);

  // tests to be done ...

  //CPPUNIT_TEST(testBuildModel);
  CPPUNIT_TEST(testBuildModel1);
  CPPUNIT_TEST(testBuildModel2);
  CPPUNIT_TEST(testSetStrategyPtr);
  CPPUNIT_TEST(testsetNonSmoothDynamicalSystemPtr);
  CPPUNIT_TEST(testsetSiconosModelXMLPtr);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildModel1();
  void testBuildModel2();
  void testSetStrategyPtr();
  void testsetNonSmoothDynamicalSystemPtr();
  void testsetSiconosModelXMLPtr();

  // Members

  double t0, T;

public:
  void setUp();
  void tearDown();

};

#endif




