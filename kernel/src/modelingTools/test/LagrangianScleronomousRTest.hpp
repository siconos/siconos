/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#ifndef __LagrangianScleronomousRTest__
#define __LagrangianScleronomousRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LagrangianScleronomousR.hpp"
#include "NonSmoothDynamicalSystem.hpp"

class LagrangianScleronomousRTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianScleronomousRTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianScleronomousRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianScleronomousR2);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianScleronomousR0();
  void testBuildLagrangianScleronomousR2();

public:
  void setUp();
  void tearDown();

};

#endif




