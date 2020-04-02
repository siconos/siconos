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
#ifndef __SiconosMemoryTest__
#define __SiconosMemoryTest__

/* /\* #include "SiconosVector.hpp" *\/ */
#include "SiconosMemory.hpp"
#include <cppunit/extensions/HelperMacros.h>

class SiconosMemoryTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosMemoryTest);


  // Test suite
  CPPUNIT_TEST_SUITE(SiconosMemoryTest);

  CPPUNIT_TEST(testBuildMemory1);
  // CPPUNIT_TEST(testBuildMemory2);
  // CPPUNIT_TEST(testBuildMemory3);
  CPPUNIT_TEST(testSetVectorMemory);
  CPPUNIT_TEST(testGetSiconosVector);
  CPPUNIT_TEST(testSwap);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_SUITE_END();

  void testBuildMemory1();
  // void testBuildMemory2();
  // void testBuildMemory3();
  void testSetVectorMemory();
  void testGetSiconosVector();
  void testSwap();
  void End();

  SP::MemoryContainer V1, V2, V3;
  SP::SiconosVector q1, q2, q3;
  SP::BlockVector c1, c2;
  unsigned int _sizeMem;
public:
  void setUp();
  void tearDown();

};

#endif
