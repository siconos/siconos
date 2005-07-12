//$id$

#ifndef PLATFORMTEST_H
#define PLATFORMTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "Model.h"
#include "LCP.h"
#include "RuntimeException.h"
#include "XMLException.h"

using namespace std;

class PlatformTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PlatformTest);

  CPPUNIT_TEST(testManualCreation1);
  CPPUNIT_TEST(testManualCreation2);
  CPPUNIT_TEST(testManualCreation3);
  CPPUNIT_TEST_EXCEPTION(testPlatformException, RuntimeException);

  CPPUNIT_TEST_SUITE_END();

  void testManualCreation1();
  void testManualCreation2();
  void testManualCreation3();
  void testPlatformException();

  //void testFail();

public:
  void setUp();
  void tearDown();
};

#endif // PLATFORMTEST_H

