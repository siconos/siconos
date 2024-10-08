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
#include "LagrangianLinearTIDSTest.hpp"

#include "BlockMatrix.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearTIDSTest);

void LagrangianLinearTIDSTest::setUp() {
  q0(0) = 1;
  q0(1) = 2;
  q0(2) = 3;
  velocity0(0) = 4;
  velocity0(1) = 5;
  velocity0(2) = 6;

  mass.setZero();
  mass(0, 0) = 1;
  mass(1, 1) = 2;
  mass(2, 2) = 3;
  K = siconos::algebra::readMatrixFromFile("K.dat");
  C = siconos::algebra::readMatrixFromFile("C.dat");

  minus_inv_M.setZero();
  minus_inv_M(0, 0) = -1.;
  minus_inv_M(1, 1) = -0.5;
  minus_inv_M(2, 2) = -1. / 3.;
  rhsK = minus_inv_M * K;
  rhsC = minus_inv_M * C;
}

void LagrangianLinearTIDSTest::tearDown() {}

// Mass, K, C
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto ds = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, velocity0, mass);
  ds->setStiffnessMatrix(K);
  ds->setDampingMatrix(C);
  siconos::algebra::SiconosVector zero(3);
  zero.setZero();

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", Type::value(*ds) ==
  // Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->q()) == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->velocity()) == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", ds->acceleration() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->mass() == mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->stiffnessMatrix() == K,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->dampingMatrix() == C,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", ds->hasExternalForces() == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->LUMass() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", ds->jacobianTotalForcesOver_q() == K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", ds->jacobianTotalForcesOver_velocity() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", ds->computeKineticEnergy() == 87.0, true);

  double time = 1.;
  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);


  
  siconos::algebra::SiconosVector acc0 = K * q0;
  acc0 += C * velocity0;
  acc0 += minus_inv_M * acc0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(rhs0, velocity0, acc0);

  auto m0 = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  m0->setZero();
  siconos::algebra::SiconosMatrix i0(3, 3);
  i0.setZero();
  i0(0, 0) = i0(1, 1) = i0(2, 2) = 1.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->n() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->x0()) == x0, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->rhs()) == rhs0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0, 1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1, 0)) == rhsK, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1, 1)) == rhsC, true);
  siconos::algebra::SiconosVector f0(3);
  f0 = K * q0 + C * velocity0;
  auto external_forces = [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result << 0, time, 2. * time;
  };

  ds->setComputeFextFunction(external_forces);

  time = 1.5;
  ds->computeTotalForces(velocity0, q0, time);
  auto x01 = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  f0 = f0 - (time * *x01);
  //   add(f0, -time * *x01, f0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->fext() == time * *x01,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", ds->totalForces() == -1. * f0, true);
  ds->computeRhs(time);
  f0 = minus_inv_M * f0;
  siconos::algebra::SiconosVector rhs1;
  siconos::algebra::concatenateVectors(rhs1, velocity0, f0);
  ds->computeJacobianRhsx(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->rhs()) == rhs1,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0, 1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1, 0)) == rhsK, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1, 1)) == rhsC, true);
  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

// Initial conditions and mass
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS2() {
  std::cout << "--> Test: constructor 2." << std::endl;
  auto ds = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, velocity0, mass);
  siconos::algebra::SiconosVector zero(3);
  zero.setZero();

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", Type::value(*ds) ==
  // Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->q()) == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", *(ds->velocity()) == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", ds->acceleration() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->mass() == mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", ds->hasExternalForces() == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->LUMass() == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", ds->computeKineticEnergy() == 87.0, true);

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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->n() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->x0()) == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->rhs()->isApprox(rhs0),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(0, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(0, 1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(1, 0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(1, 1)) == *m0, true);
  siconos::algebra::SiconosVector f0(3);

  auto external_forces = [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
    auto count = 0;
    for (auto& v : result) v = count++ * time;
  };

  ds->setComputeFextFunction(external_forces);

  time = 1.5;
  ds->computeTotalForces(velocity0, q0, time);
  auto x01 = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  //   add(f0, -time * *x01, f0);
  f0 = f0 - (time * *x01);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->fext() == time * *x01,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", ds->totalForces().isApprox(-1. * f0), true);
  std::cout << "--> Constructor 2 test ended with success." << std::endl;
}
