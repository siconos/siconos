/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef MultiBodyTest_h
#define MultiBodyTest_h

#include <cppunit/extensions/HelperMacros.h>

class MultiBodyTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(MultiBodyTest);

  // tests to be done ...
  CPPUNIT_TEST(t1);

  CPPUNIT_TEST(t2);

  //  CPPUNIT_TEST(t3);

  //  CPPUNIT_TEST(t4);

  //  CPPUNIT_TEST(t5);

  //  CPPUNIT_TEST(t6);

  //  CPPUNIT_TEST(t7);
  //  CPPUNIT_TEST(t8);

  CPPUNIT_TEST_SUITE_END();

  // Members
  void t1();
  void t2();
  void t3();
  void t4();
  void t5();
  void t6();
  void t7();
  void t8();

public:
  void setUp();
  void tearDown();

};

#endif
