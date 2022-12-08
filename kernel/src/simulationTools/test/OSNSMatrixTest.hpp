/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef __OSNSMatrixTest__
#define __OSNSMatrixTest__

#include <cppunit/extensions/HelperMacros.h>
#include "OSNSMatrix.hpp"
#include "TimeStepping.hpp"
#include "OneStepNSProblem.hpp"

class OSNSMatrixTest : public CppUnit::TestFixture
{

private:
  
  ACCEPT_SERIALIZATION(OSNSMatrixTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(OSNSMatrixTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildOSNSMatrix0);
  CPPUNIT_TEST(testBuildOSNSMatrix1);
  CPPUNIT_TEST(testBuildOSNSMatrix2);
  CPPUNIT_TEST(testFill);
  CPPUNIT_TEST(testConvert);
  CPPUNIT_TEST(testFill2);
  CPPUNIT_TEST(testConvert2);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildOSNSMatrix0();
  void testBuildOSNSMatrix1();
  void testBuildOSNSMatrix2();
  void testFill();
  void testConvert();
  void testFill2();
  void testConvert2();
  void End();
  // Members

  unsigned int n;
  double tol;
  SP::InteractionsGraph indexSet;
  MapOfMapOfInteractionMatrices blocks;
  Model * temp ;
public:
  void setUp();
  void tearDown();

};

#endif




