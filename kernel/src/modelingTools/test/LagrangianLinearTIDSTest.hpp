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
#ifndef __LagrangianLinearTIDSTest__
#define __LagrangianLinearTIDSTest__

#include <cppunit/extensions/HelperMacros.h>
#include "LagrangianLinearTIDS.hpp"
#include "RuntimeException.hpp"

class LagrangianLinearTIDSTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDSTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianLinearTIDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianLinearTIDS1);
  CPPUNIT_TEST(testBuildLagrangianLinearTIDS2);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianLinearTIDS1();
  void testBuildLagrangianLinearTIDS2();

  // Members

  SP::SiconosVector q0, velocity0;
  SP::SiconosMatrix mass, K, C, rhsK, rhsC, minus_inv_M;

public:
  void setUp();
  void tearDown();

};

#endif




