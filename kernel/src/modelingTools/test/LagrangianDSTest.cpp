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
#include "LagrangianDSTest.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
#include "BlockMatrix.hpp"
#include "SiconosMatrixOp.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianDSTest);

void LagrangianDSTest::setUp()
{
  q0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*q0)(0) = 1;
  (*q0)(1) = 2;
  (*q0)(2) = 3;
  velocity0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*velocity0)(0) = 4;
  (*velocity0)(1) = 5;
  (*velocity0)(2) = 6;

  u0 = std::make_shared<siconos::algebra::SiconosVector>(2);
  (*u0)(0) = 4;
  (*u0)(1) = 5;

  mass = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  mass->setZero();
  (*mass)(0, 0) = 1;
  (*mass)(1, 1) = 2;
  (*mass)(2, 2) = 3;
}

void LagrangianDSTest::tearDown() {}

// constructor from initial state only
void LagrangianDSTest::testBuildLagrangianDS1()
{
  std::cout << "--> Test: constructor 1." << std::endl;

  siconos::algebra::SiconosVector zero(3);
  zero.setZero();
  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildLagrangianDS1 : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->velocity0() == velocity0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->mass() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->forces() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->fInt() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->fExt() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->fGyr() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->inverseMass() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianFIntq() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianFIntqDot() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianFGyrq() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianFGyrqDot() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianqForces() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianvForces() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->computeKineticEnergy() == 38.5,
                               true);

  double time = 1.;
  auto m0 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  m0->setZero();

  ds->computeForces(time, ds->q(), ds->velocity());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->forces()->isApprox(zero), true);
  ds->computeJacobianqForces(time);
  ds->computeJacobianqDotForces(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianqForces() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->jacobianvForces() == nullptr,
                               true);

  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, *q0, *velocity0);
  siconos::algebra::concatenateVectors(rhs0, *velocity0, zero);
  siconos::algebra::SiconosMatrix i0(3, 3);
  i0(0, 0) = i0(1, 1) = i0(2, 2) = 1.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->n() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", *(ds->x0()) == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", *(ds->rhs()) == rhs0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS1 : ", *(ds->jacobianRhsx()->block(0, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS1 : ", *(ds->jacobianRhsx()->block(0, 1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS1 : ", *(ds->jacobianRhsx()->block(1, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS1 : ", *(ds->jacobianRhsx()->block(1, 1)) == *m0, true);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

// constructor from initial state and mass matrix
void LagrangianDSTest::testBuildLagrangianDS4()
{
  std::cout << "--> Test: constructor 4." << std::endl;

  siconos::algebra::SiconosVector zero(3);
  zero.setZero();
  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0, mass);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildLagrangianDS4 : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->velocity0() == velocity0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", *(ds->mass()) == *(mass), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->forces() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->fInt() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->fExt() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->fGyr() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->inverseMass() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->jacobianFIntq() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->jacobianFIntqDot() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->jacobianFGyrq() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->jacobianFGyrqDot() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->jacobianqForces() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->jacobianvForces() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildLagrangianDS : ", ds->computeKineticEnergy() == 87.0,
                               true);

  double time = 1.;
  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, *q0, *velocity0);
  siconos::algebra::concatenateVectors(rhs0, *velocity0, zero);
  auto m0 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  m0->setZero();
  siconos::algebra::SiconosMatrix i0(3, 3);
  i0.setZero();
  i0(0, 0) = i0(1, 1) = i0(2, 2) = 1.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->n() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", *(ds->x0()) == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", *(ds->rhs()) == rhs0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS4 : ", *(ds->jacobianRhsx()->block(0, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS4 : ", *(ds->jacobianRhsx()->block(0, 1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS4 : ", *(ds->jacobianRhsx()->block(1, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS4 : ", *(ds->jacobianRhsx()->block(1, 1)) == *m0, true);

  std::cout << "--> Constructor 4 test ended with success." << std::endl;
}

// constructor from initial state and plugged mass
void LagrangianDSTest::testBuildLagrangianDS5()
{
  std::cout << "--> Test: constructor 5." << std::endl;
  std::string plugin = "TestPlugin:computeMass";
  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0, plugin);

  siconos::algebra::SiconosVector zero(3);
  zero.setZero();
  auto m0 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  m0->setZero();

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildLagrangianDS5 : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->velocity0() == velocity0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", *(ds->mass()) == *m0, true);
  ds->computeMass();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", *(ds->mass()) == *(mass), true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->forces() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->fInt() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->fExt() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->fGyr() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->inverseMass() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->jacobianFIntq() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->jacobianFIntqDot() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->jacobianFGyrq() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->jacobianFGyrqDot() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->jacobianqForces() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->jacobianvForces() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->computeKineticEnergy() == 87.0,
                               true);

  double time = 1.;
  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, *q0, *velocity0);
  siconos::algebra::concatenateVectors(rhs0, *velocity0, zero);
  siconos::algebra::SiconosMatrix i0(3, 3);
  i0(0, 0) = i0(1, 1) = i0(2, 2) = 1.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->n() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", *(ds->x0()) == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", *(ds->rhs()) == rhs0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS5 : ", *(ds->jacobianRhsx()->block(0, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS5 : ", *(ds->jacobianRhsx()->block(0, 1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS5 : ", *(ds->jacobianRhsx()->block(1, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianDS5 : ", *(ds->jacobianRhsx()->block(1, 1)) == *m0, true);

  std::cout << "--> Constructor 5 test ended with success." << std::endl;
}
