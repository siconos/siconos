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
#ifndef __FirstOrderNonLinearDSTest__
#define __FirstOrderNonLinearDSTest__

#include <cppunit/extensions/HelperMacros.h>
#include "FirstOrderNonLinearDS.hpp"
#include "RuntimeException.hpp"

class FirstOrderNonLinearDSTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderNonLinearDSTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderNonLinearDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS1);
  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS2);
  CPPUNIT_TEST(testBuildFirstOrderNonLinearDS3);
  CPPUNIT_TEST(testSetX0);
  CPPUNIT_TEST(testSetX0Ptr);
  CPPUNIT_TEST(testSetx);
  CPPUNIT_TEST(testSetxPtr);
  CPPUNIT_TEST(testSetR);
  CPPUNIT_TEST(testSetRPtr);
  CPPUNIT_TEST(testSetJacobianfxPtr);
  CPPUNIT_TEST(testInitMemory);
  CPPUNIT_TEST(testSwap);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderNonLinearDS1();
  void testBuildFirstOrderNonLinearDS2();
  void testBuildFirstOrderNonLinearDS3();
  void testSetX0();
  void testSetX0Ptr();
  void testSetx();
  void testSetxPtr();
  void testSetR();
  void testSetRPtr();
  void testSetJacobianfxPtr();
  void testInitMemory();
  void testSwap();

  // Members

  SP::SiconosVector x0, xnull;
  SP::SiconosMatrix J0, M;
public:
  void setUp();
  void tearDown();

};

#endif




