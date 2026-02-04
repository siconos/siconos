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

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
#include "io.hpp"

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
  K = siconos::algebra::io::readDenseMatrix("K.dat");
  C = siconos::algebra::io::readDenseMatrix("C.dat");

  minus_inv_M.setZero();
  minus_inv_M(0, 0) = -1.;
  minus_inv_M(1, 1) = -0.5;
  minus_inv_M(2, 2) = -1. / 3.;
  rhsK = minus_inv_M * K;
  rhsC = minus_inv_M * C;
}

void LagrangianLinearTIDSTest::tearDown() {}

// Mass, K, C, alias mode
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS_alias() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto ds = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
      q0, velocity0, mass, siconos::algebra::alias_t);

  ds->setStiffnessMatrix(K, siconos::algebra::alias_t);
  ds->setDampingMatrix(C, siconos::algebra::alias_t);
  siconos::algebra::SiconosVector zero = siconos::algebra::SiconosVector::Zero(3);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", Type::value(*ds) ==
  // Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->q_read() == q0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->velocity_read() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->acceleration() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->mass() == mass,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->stiffnessMatrix() == K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->dampingMatrix() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", *(ds->p(1)) == zero,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->p(0) == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->p(2) == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->hasExternalForces() == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->LUMass() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->computeKineticEnergy() == 87.0, true);

  // Check alias (shared memory)
  K *= 2;
  C *= 3;
  q0 *= 2;
  mass *= 6.2;
  velocity0 *= 3;
  ds->resetToInitialState();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->mass() == mass,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->stiffnessMatrix() == K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->dampingMatrix() == C, true);

  double time = 1.;
  ds->initRhs(time);

  siconos::algebra::SiconosVector x0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);
  minus_inv_M = -mass.inverse();
  siconos::algebra::SiconosVector acc0 = minus_inv_M * (K * q0 + C * velocity0);
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(rhs0, velocity0, acc0);

  siconos::algebra::SiconosMatrix m0 = siconos::algebra::SiconosMatrix::Zero(3, 3);
  siconos::algebra::SiconosMatrix i0(3, 3);
  i0.setZero();
  i0(0, 0) = i0(1, 1) = i0(2, 2) = 1.;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->x_size() == 2 * 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", ds->x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->rhs()->isApprox(rhs0), true);

  siconos::algebra::SiconosVector jacoref{36};
  jacoref.setZero();

  for (unsigned int j = 0; j < 3; ++j) {
    jacoref((3 + j) * 6 + j) = 1.0;
  }
  rhsK = minus_inv_M * K;
  rhsC = minus_inv_M * C;

  for (unsigned int j = 0; j < 3; ++j) {
    for (unsigned int i = 0; i < 3; ++i) jacoref(j * 6 + i + 3) = rhsK(i, j);
    for (unsigned int i = 0; i < 3; ++i) jacoref((j + 3) * 6 + i + 3) = rhsC(i, j);
  }
  const auto& jaco_rhs = ds->jacobianRhsOver_x();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 11: ", jacoref.isApprox(jaco_rhs, 1e-12),
                               true);
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
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->fext() == time * *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_alias : ", ds->totalForces() == -1. * f0, true);
  ds->computeRhs(time);
  f0 = minus_inv_M * f0;
  siconos::algebra::SiconosVector rhs1;
  siconos::algebra::concatenateVectors(rhs1, velocity0, f0);
  ds->computeJacobianRhsOver_x(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_alias : ", *(ds->rhs()) == rhs1,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 11: ", jacoref.isApprox(jaco_rhs, 1e-12),
                               true);
  std::cout << "✅ Basic constructor (alias) test ended with success." << std::endl;
}

void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS_copy() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto ds = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
      q0, velocity0, mass, siconos::algebra::copy_t);

  ds->setStiffnessMatrix(K, siconos::algebra::copy_t);
  ds->setDampingMatrix(C, siconos::algebra::copy_t);
  siconos::algebra::SiconosVector zero = siconos::algebra::SiconosVector::Zero(3);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", Type::value(*ds) ==
  // Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->q_read() == q0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->velocity_read() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->acceleration() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->mass() == mass,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->stiffnessMatrix() == K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->dampingMatrix() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", *(ds->p(1)) == zero,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->p(0) == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->p(2) == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->hasExternalForces() == false, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->LUMass() == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->computeKineticEnergy() == 87.0, true);

  // Check alias (shared memory)
  K *= 2;
  C *= 3;
  q0 *= 2;
  mass *= 6.2;
  velocity0 *= 3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->q0() == q0, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->velocity0() == velocity0, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS_copy : ", ds->mass() == mass,
                               false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->stiffnessMatrix() == K, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS_copy : ", ds->dampingMatrix() == C, false);
  std::cout << "✅ Basic constructor (copy) test ended with success." << std::endl;
}

// Initial conditions and mass
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS2() {
  std::cout << "--> Test: constructor 2." << std::endl;
  auto ds = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
      q0, velocity0, mass, siconos::algebra::alias_t);
  siconos::algebra::SiconosVector zero = siconos::algebra::SiconosVector::Zero(3);

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

  siconos::algebra::SiconosMatrix m0 = siconos::algebra::SiconosMatrix::Zero(3, 3);
  siconos::algebra::SiconosMatrix i0(3, 3);
  i0.setZero();
  i0(0, 0) = i0(1, 1) = i0(2, 2) = 1.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->x_size() == 2 * 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->rhs()->isApprox(rhs0),
                               true);
  siconos::algebra::SiconosVector jacoref{36};
  jacoref.setZero();
  for (unsigned int j = 0; j < 3; ++j) {
    jacoref((3 + j) * 6 + j) = 1.0;
  }
  const auto& jaco_rhs = ds->jacobianRhsOver_x();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 11: ", jacoref.isApprox(jaco_rhs, 1e-12),
                               true);
  auto external_forces = [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
    auto count = 0;
    for (auto& v : result) v = count++ * time;
  };

  ds->setComputeFextFunction(external_forces);

  time = 1.5;

  ds->computeTotalForces(velocity0, q0, time);

  siconos::algebra::SiconosVector3 ref;
  ref << 0, 1, 2;
  ref *= time;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->fext() == ref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianLinearTIDS2 : ", ds->totalForces().isApprox(ref), true);
  std::cout << "✅ constructor test (2) ended with success." << std::endl;
}
