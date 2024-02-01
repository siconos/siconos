/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef KERNEL_TEST_HPP
#define KERNEL_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "SiconosConfig.h"

class KernelTest : public CppUnit::TestFixture {
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

  CPPUNIT_TEST(t9);

  CPPUNIT_TEST_SUITE_END();

  void t0();
  void t1();
  void t2();
  void t3();
  // void t4();
  void t5();
  void t6();

#ifdef HAVE_SICONOS_MECHANICS
  void t7();
  void t8();
#endif

  void t9();

  std::string BBxml;

 public:
  void setUp();
  void tearDown();
};

#endif
