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
#include "NewtonEulerDSTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NewtonEulerDSTest);

void NewtonEulerDSTest::setUp() {
  q0 << 1, 2, 3, 0., 1, 0, 0;
  q01 << 1, 2, 3, 1, 0, 0, 0.;
  twist0 << 4., 5, 6, 7, 8, 9;
  inertia(0, 0) = 1;
  inertia(1, 1) = 2;
  inertia(2, 2) = 3;
}

void NewtonEulerDSTest::tearDown() {}

// constructor from data
void NewtonEulerDSTest::testBuildNewtonEulerDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;

  siconos::modeling::NewtonEulerDS neds{q0, twist0, mass, inertia};

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildNewtonEulerDS1A : ", Type::value(*ds) == Type::NewtonEulerDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1B : ", neds.number() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", neds.dimension() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", neds.getqDim() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", neds.scalarMass() == mass, true);

  siconos::algebra::SiconosMatrix massMatrix{6, 6};
  massMatrix(0, 0) = mass;
  massMatrix(1, 1) = mass;
  massMatrix(2, 2) = mass;
  massMatrix.block<3, 3>(3, 3) = inertia;

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildNewtonEulerDS1D : ", neds.totalInertiaMatrix() == massMatrix, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildNewtonEulerDS1D : ", neds.computeKineticEnergy() == 595.0, true);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}
// constructor from data
