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

#include "BlockMatrix.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianDSTest);

void LagrangianDSTest::setUp() {
  q0 << 1, 2, 3;
  velocity0 << 4, 5, 6;

  mass.setZero();
  mass(0, 0) = 1;
  mass(1, 1) = 2;
  mass(2, 2) = 3;
}

void LagrangianDSTest::tearDown() {}

// constructor from initial state only
void LagrangianDSTest::testBuildLagrangianDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  siconos::algebra::SiconosVector zero = siconos::algebra::SiconosVector::Zero(3);
  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildLagrangianDS1 : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->velocity0() == velocity0,
                               true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS1 : ", ds->computeKineticEnergy() == 38.5,
                               true);

  double time = 1.;
  auto m0 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  m0->setZero();

  // Try to compute forces and jacobians, even if not set, just to check
  // that nothing happens and that no exception is raised
  ds->computeTotalForces(ds->velocity_read(), ds->q_read(), time);
  ds->computeJacobianTotalForcesOver_q(ds->velocity_read(), ds->q_read(), time);
  ds->computeJacobianTotalForcesOver_velocity(ds->velocity_read(), ds->q_read(), time);

  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);
  siconos::algebra::concatenateVectors(rhs0, velocity0, zero);
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
void LagrangianDSTest::testBuildLagrangianDS4() {
  std::cout << "--> Test: constructor 4." << std::endl;

  siconos::algebra::SiconosVector zero(3);
  zero.setZero();
  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0);
  ds->setConstantMass(mass);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildLagrangianDS4 : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->velocity0() == velocity0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->mass() == mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->p_read(1) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildLagrangianDS : ", ds->computeKineticEnergy() == 87.0,
                               true);

  double time = 1.;
  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);
  siconos::algebra::concatenateVectors(rhs0, velocity0, zero);
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
void LagrangianDSTest::testBuildLagrangianDS5() {
  std::cout << "--> Test: constructor 5." << std::endl;

  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0);

  //   auto mass_func = [](const Eigen::Ref<const siconos::algebra::SiconosVector>& pos,
  //                       Eigen::Ref<siconos::algebra::MapType> result) {
  auto mass_func = [](const Eigen::Ref<const siconos::algebra::SiconosVector>& pos,
                      Eigen::Ref<siconos::algebra::MapType> result) {
    result.setZero();
    result(0, 0) = 0. * pos(0);  // just to check that pos is callable ...
    result(0, 0) = 1;
    result(1, 1) = 2.;
    result(2, 2) = 3.;
  };

  ds->setComputeMassFunction(mass_func);
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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->mass() == *m0, true);
  ds->computeMass(*ds->q());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->mass() == mass, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5 : ", ds->computeKineticEnergy() == 87.0,
                               true);

  double time = 1.;
  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);
  siconos::algebra::concatenateVectors(rhs0, velocity0, zero);
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
