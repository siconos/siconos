
#pragma once

#include <cppunit/extensions/HelperMacros.h>

class RtOsiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(RtOsiTest);
  CPPUNIT_TEST(testOsi0);
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp() override;
  void tearDown() override;
  void testOsi0();
};
