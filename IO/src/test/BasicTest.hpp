#ifndef BASIC_TEST_HPP
#define BASIC_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>

class BasicTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(BasicTest);

  CPPUNIT_TEST(t0);
  CPPUNIT_TEST(t1);
  CPPUNIT_TEST(t2);
  CPPUNIT_TEST(t3);
  CPPUNIT_TEST(t4);

  CPPUNIT_TEST_SUITE_END();

  void t0();
  void t1();
  void t2();
  void t3();
  void t4();

public:
  void setUp();
  void tearDown();
};

#endif
