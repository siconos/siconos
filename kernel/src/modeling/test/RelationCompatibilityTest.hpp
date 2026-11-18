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
#ifndef __RelationCompatibilityTest__
#define __RelationCompatibilityTest__

#include <cppunit/extensions/HelperMacros.h>

/** Checks the deprecated names of the LagrangianR / NewtonEulerR API (kept for one release):
 *  - a class overriding a former virtual name (e.g. computeJacobianhOver_q) is still called
 *    through the new name (computeH), which is the one used by the kernel;
 *  - a class overriding the new name is reached through the former name;
 *  - without override, both names run the default implementation (user function);
 *  - the deprecated accessors, setters, members and class names forward to the new ones.
 */
class RelationCompatibilityTest : public CppUnit::TestFixture {
 private:
  CPPUNIT_TEST_SUITE(RelationCompatibilityTest);

  CPPUNIT_TEST(testLagrangianScleronomousR);
  CPPUNIT_TEST(testLagrangianRheonomousR);
  CPPUNIT_TEST(testLagrangianCompliantR);
  CPPUNIT_TEST(testNewtonEulerR);
  CPPUNIT_TEST(testClassAliases);
  CPPUNIT_TEST_SUITE_END();

  void testLagrangianScleronomousR();
  void testLagrangianRheonomousR();
  void testLagrangianCompliantR();
  void testNewtonEulerR();
  void testClassAliases();

 public:
  void setUp() {}
  void tearDown() {}
};

#endif
