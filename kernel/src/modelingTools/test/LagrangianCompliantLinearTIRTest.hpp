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
#ifndef __LagrangianCompliantLinearTIRTest__
#define __LagrangianCompliantLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "NonSmoothDynamicalSystem.hpp"
#include "LagrangianCompliantLinearTIR.hpp"

class LagrangianCompliantLinearTIRTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianCompliantLinearTIRTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(LagrangianCompliantLinearTIRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildLagrangianCompliantLinearTIR1);
  CPPUNIT_TEST(testBuildLagrangianCompliantLinearTIR2);
  CPPUNIT_TEST(testBuildLagrangianCompliantLinearTIR3);
  CPPUNIT_TEST(testSetCPtr);
  CPPUNIT_TEST(testSetDPtr);
  CPPUNIT_TEST(testSetFPtr);
  CPPUNIT_TEST(testSetEPtr);
  CPPUNIT_TEST(testGetJacPtr);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildLagrangianCompliantLinearTIR0();
  void testBuildLagrangianCompliantLinearTIR1();
  void testBuildLagrangianCompliantLinearTIR2();
  void testBuildLagrangianCompliantLinearTIR3();

  void testSetCPtr();
  void testSetDPtr();
  void testSetFPtr();
  void testSetEPtr();
  void testGetJacPtr();
  void End();

  // Members

  SP::SimpleMatrix C, B, F, D;
  SP::SiconosVector e;
  SP::NonSmoothDynamicalSystem nsds;

public:
  void setUp();
  void tearDown();

};

#endif




