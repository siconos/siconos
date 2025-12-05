
#pragma once

#include <cppunit/extensions/HelperMacros.h>

class StorageTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(StorageTest);
  CPPUNIT_TEST(testStorage0);
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp() override;
  void tearDown() override;
  void testStorage0();
};
