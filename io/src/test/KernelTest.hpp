#ifndef KERNEL_TEST_HPP
#define KERNEL_TEST_HPP

#include "SiconosConfig.h"
#include <cppunit/extensions/HelperMacros.h>

class KernelTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(KernelTest);

  CPPUNIT_TEST(t0);
  CPPUNIT_TEST(t1);
  CPPUNIT_TEST(t2);
  CPPUNIT_TEST(t3);
//  CPPUNIT_TEST(t4);
  CPPUNIT_TEST(t5);
  CPPUNIT_TEST(t6);

#ifdef HAVE_SICONOS_MECHANICS
  CPPUNIT_TEST(t7);
  CPPUNIT_TEST(t8);
#endif

  CPPUNIT_TEST_SUITE_END();

  void t0();
  void t1();
  void t2();
  void t3();
  //void t4();
  void t5();
  void t6();

#ifdef HAVE_SICONOS_MECHANICS
  void t7();
  void t8();
#endif

  std::string BBxml;
public:
  void setUp();
  void tearDown();
};

#endif
