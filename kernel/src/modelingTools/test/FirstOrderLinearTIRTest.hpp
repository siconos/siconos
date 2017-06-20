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
#ifndef __FirstOrderLinearTIRTest__
#define __FirstOrderLinearTIRTest__

#include <cppunit/extensions/HelperMacros.h>
#include "NonSmoothDynamicalSystem.hpp"
#include "FirstOrderLinearTIR.hpp"

class FirstOrderLinearTIRTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(FirstOrderLinearTIRTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearTIRTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearTIR1);
  CPPUNIT_TEST(testBuildFirstOrderLinearTIR2);
  //  CPPUNIT_TEST(testSetC);
  CPPUNIT_TEST(testSetCPtr);
  //  CPPUNIT_TEST(testSetD);
  CPPUNIT_TEST(testSetDPtr);
  //  CPPUNIT_TEST(testSetF);
  CPPUNIT_TEST(testSetFPtr);
  //  CPPUNIT_TEST(testSetE);
  CPPUNIT_TEST(testSetEPtr);
  //  CPPUNIT_TEST(testSetB);
  CPPUNIT_TEST(testSetBPtr);
  //  CPPUNIT_TEST(testGetJac);
  CPPUNIT_TEST(testGetJacPtr);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearTIR0();
  void testBuildFirstOrderLinearTIR1();
  void testBuildFirstOrderLinearTIR2();
  //  void testSetC();
  void testSetCPtr();
  //  void testSetD();
  void testSetDPtr();
  //  void testSetF();
  void testSetFPtr();
  //  void testSetE();
  void testSetEPtr();
  //  void testSetB();
  void testSetBPtr();
  //  void testGetJac();
  void testGetJacPtr();

  // Members

  SP::SimpleMatrix C, B, F, D;
  SP::SiconosVector e;
  SP::NonSmoothDynamicalSystem nsds;

public:
  void setUp();
  void tearDown();

};

#endif




