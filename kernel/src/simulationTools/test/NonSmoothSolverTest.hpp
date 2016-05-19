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
#ifndef __NonSmoothSolverTest__
#define __NonSmoothSolverTest__

#include <cppunit/extensions/HelperMacros.h>
#include "NonSmoothSolver.hpp"
#include "LagrangianLinearTIDS.hpp"

class NonSmoothSolverTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NonSmoothSolverTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(NonSmoothSolverTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildNonSmoothSolver0);
  CPPUNIT_TEST(testBuildNonSmoothSolver1);
  CPPUNIT_TEST(testBuildNonSmoothSolver2);
  CPPUNIT_TEST(testBuildNonSmoothSolver3);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildNonSmoothSolver0();
  void testBuildNonSmoothSolver1();
  void testBuildNonSmoothSolver2();
  void testBuildNonSmoothSolver3();
  void End();
  // Members

  std::vector<int> iparam;
  std::vector<double> dparam;

public:
  void setUp();
  void tearDown();

};

#endif




