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
#ifndef __FirstOrderLinearDSTest__
#define __FirstOrderLinearDSTest__

#include <cppunit/extensions/HelperMacros.h>
#include "FirstOrderLinearDS.hpp"
#include "RuntimeException.hpp"

class FirstOrderLinearDSTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearDSTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(FirstOrderLinearDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildFirstOrderLinearDS0);
  CPPUNIT_TEST(testBuildFirstOrderLinearDS1);
  //  CPPUNIT_TEST(testSetA);
  CPPUNIT_TEST(testSetAPtr);
  //  CPPUNIT_TEST(testSetB);
  CPPUNIT_TEST(testSetBPtr);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildFirstOrderLinearDS0();
  void testBuildFirstOrderLinearDS1();
  //  void testSetA();
  void testSetAPtr();
  //  void testSetB();
  void testSetBPtr();
  // Members

  SP::SiconosVector x0;
  SP::SiconosVector b0;
  SP::SiconosMatrix A0;
public:
  void setUp();
  void tearDown();

};

#endif




