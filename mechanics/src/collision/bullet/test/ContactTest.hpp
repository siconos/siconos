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
#ifndef TestContact_h
#define TestContact_h

#include <cppunit/extensions/HelperMacros.h>

class ContactTest : public CppUnit::TestFixture
{
private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(ContactTest);

  // tests to be done ...
  CPPUNIT_TEST(t1);
  CPPUNIT_TEST(t2);
  CPPUNIT_TEST(t3);
  CPPUNIT_TEST(t4);

  CPPUNIT_TEST_SUITE_END();

  // Members
  void t1();
  void t2();
  void t3();
  void t4();

public:
  void setUp();
  void tearDown();
};

#endif
