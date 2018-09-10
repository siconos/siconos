/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef __EigenProblemsTest__
#define __EigenProblemsTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SimpleMatrix.hpp"
#include "EigenProblems.hpp"

class EigenProblemsTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(EigenProblemsTest);

  // test suite
  CPPUNIT_TEST_SUITE(EigenProblemsTest);

  CPPUNIT_TEST(testSyev);
  CPPUNIT_TEST(testGeev1);
  CPPUNIT_TEST(testGeev2);
  CPPUNIT_TEST(testGeev3);

  CPPUNIT_TEST_SUITE_END();

  void testSyev();
  void testGeev1();
  void testGeev2();
  void testGeev3();

  void End();

  unsigned int size;
  SP::SiconosMatrix A, Aref;

public:
  void setUp();
  void tearDown();

};

#endif
