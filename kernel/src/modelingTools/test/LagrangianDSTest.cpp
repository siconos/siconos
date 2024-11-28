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
  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildLagrangianDS1: ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 1: ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 2: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 3: ", ds->velocity0() == velocity0, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 4: ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 5: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 6: ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 7: ", ds->computeKineticEnergy() == 38.5,
                               true);

  double time = 1.;

  // Try to compute forces and jacobians, even if not set, just to check
  // that nothing happens and that no exception is raised
  ds->computeTotalForces(ds->velocity_read(), ds->q_read(), time);
  ds->computeJacobianTotalForcesOver_q(ds->velocity_read(), ds->q_read(), time);
  ds->computeJacobianTotalForcesOver_velocity(ds->velocity_read(), ds->q_read(), time);

  // Init and compute RHS (First-order formulation)
  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);
  siconos::algebra::SiconosVector zero = siconos::algebra::SiconosVector::Zero(3);
  siconos::algebra::concatenateVectors(rhs0, velocity0, zero);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 8: ", ds->x_size() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 9: ", ds->x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 1 - 10: ", *(ds->rhs()) == rhs0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 1 - 11: ", ds->jacobianRhsOver_x()->block(0, 0)->isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 1 - 12: ", ds->jacobianRhsOver_x()->block(0, 1)->isIdentity(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 1 - 13: ", ds->jacobianRhsOver_x()->block(1, 0)->isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 1 - 14: ", ds->jacobianRhsOver_x()->block(1, 1)->isZero(), true);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

// A DS with all operators (mass, fext, fint, fgyr)
void LagrangianDSTest::testBuildLagrangianDS2() {
  std::cout << "--> Test: constructor 2." << std::endl;

  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0);
  ds->setConstantMass(mass);

  siconos::algebra::SiconosVector3 fext;
  fext << 1, 2, 3;
  ds->setConstantFext(fext);

  //   siconos::algebra::SiconosVector3 fgyr;
  //   fgyr << 1, 2, 3;
  //   ds->setConstantFgyr(fgyr);

  siconos::algebra::SiconosVector3 ref;
  ref << 1, 2, 3;

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "test - Constr 2: ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->mass() == mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->computeKineticEnergy() == 87.0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->fext().isApprox(ref), true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildLagrangianDS: ", ds->fgyr().isApprox(ref),
  //  true);

  ds->setComputeFintFunction([](const siconos::algebra::SiconosVector &velocity,
                                const siconos::algebra::SiconosVector &q, double time,
                                Eigen::Ref<siconos::algebra::MapVectorType> result) {
    auto ndof = result.size();
    for (size_t i = 0; i < ndof; ++i) result[i] = velocity(i) * time + q(i);
  });

  siconos::algebra::SiconosVector3 q;
  q << 1, 2, 3;
  siconos::algebra::SiconosVector3 vel;
  vel << 4, 5, 6;
  double time = 2.;
  ds->computeFint(vel, q, time);

  siconos::algebra::SiconosMatrix33 jacq, jacvel;
  // set 'random' value for jacobians, whatever fint is, just for tests
  jacq.setConstant(1.);
  jacvel.setConstant(2.);
  ds->setConstantJacobianFintOver_q(jacq);
  ds->setConstantJacobianFintOver_velocity(jacvel);

  siconos::algebra::SiconosVector3 reffint;
  reffint << 9, 12, 15;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->fint().isApprox(reffint), true);

  ds->computeTotalForces(vel, q, time);
  ds->computeJacobianTotalForcesOver_q(vel, q, time);
  ds->computeJacobianTotalForcesOver_velocity(vel, q, time);

  siconos::algebra::SiconosMatrix33 refjac;
  refjac.setConstant(-1);
  siconos::algebra::SiconosVector3 refforces;
  refforces = fext - reffint;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->totalForces().isApprox(refforces),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 2: ", ds->jacobianTotalForcesOver_q().isApprox(refjac), true);
  refjac.setConstant(-2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 2: ", ds->jacobianTotalForcesOver_velocity().isApprox(refjac), true);

  ds->initRhs(time);
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);
  siconos::algebra::SiconosVector3 ref0;
  ref0 = mass.inverse() * refforces;
  siconos::algebra::concatenateVectors(rhs0, velocity0, ref0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->x_size() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", *(ds->rhs()) == rhs0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 2: ", ds->jacobianRhsOver_x()->block(0, 0)->isZero(),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 2: ", ds->jacobianRhsOver_x()->block(0, 1)->isIdentity(), true);
  siconos::algebra::SiconosMatrix33 refj2;
  refjac.setConstant(-1);

  refj2 = mass.inverse() * refjac;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 2: ", ds->jacobianRhsOver_x()->block(1, 0)->isApprox(refj2), true);
  refjac.setConstant(-2);
  refj2 = mass.inverse() * refjac;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 2: ", ds->jacobianRhsOver_x()->block(1, 1)->isApprox(refj2), true);

  std::cout << "--> Constructor 2 test ended with success." << std::endl;
}

// A DS with all operators (mass, fext, fint, fgyr and jacobians) as user-defined functions
void LagrangianDSTest::testBuildLagrangianDS3() {
  std::cout << "--> Test: constructor 3." << std::endl;

  auto ds = std::make_shared<siconos::modeling::LagrangianDS>(q0, velocity0);

  auto mass_func = [](const Eigen::Ref<const siconos::algebra::SiconosVector> &pos,
                      Eigen::Ref<siconos::algebra::MapType> result) {
    result.setZero();
    result(0, 0) = 0. * pos(0);  // just to check that pos is callable ...
    result(0, 0) = 1;
    result(1, 1) = 2.;
    result(2, 2) = 3.;
  };

  ds->setComputeMassFunction(mass_func);

  ds->setComputeFextFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto i = 0;
        for (auto &v : result) v = time;
      });

  ds->setComputeFgyrFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &vel,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
         Eigen::Ref<siconos::algebra::MapVectorType> result) { result = vel + q; });

  // set 'random' value for jacobians, whatever fgyr is, just for tests
  ds->setComputeJacobianFgyrOver_qFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
         Eigen::Ref<siconos::algebra::MapType> result) { result.setConstant(2.); });

  ds->setComputeJacobianFgyrOver_velocityFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
         Eigen::Ref<siconos::algebra::MapType> result) { result.setConstant(3.); });

  ds->setComputeFintFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
                                const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
                                double time,
                                Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result = v + q;
    result *= time;
  });

  ds->setComputeJacobianFintOver_qFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapType> result) { result.setConstant(4.); });

  ds->setComputeJacobianFintOver_velocityFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapType> result) { result.setConstant(5.); });

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->p(2) == nullptr, true);

  siconos::algebra::SiconosVector3 q;
  q << 1, 2, 3;
  siconos::algebra::SiconosVector3 vel;
  vel << 4, 5, 6;
  double time = 2.;

  ds->computeMass(q);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->mass() == mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->computeKineticEnergy() == 87.0, true);

  ds->computeTotalForces(vel, q, time);
  ds->computeJacobianTotalForcesOver_q(vel, q, time);
  ds->computeJacobianTotalForcesOver_velocity(vel, q, time);

  siconos::algebra::SiconosVector3 refforces;
  refforces.setConstant(time);
  refforces -= (time + 1) * (vel + q);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Constr 3: ", ds->totalForces().isApprox(refforces),
                               true);

  siconos::algebra::SiconosMatrix33 jacq, jacvel;

  siconos::algebra::SiconosMatrix33 refjac;
  refjac.setConstant(-6);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 2: ", ds->jacobianTotalForcesOver_q().isApprox(refjac), true);
  refjac.setConstant(-8);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Constr 2: ", ds->jacobianTotalForcesOver_velocity().isApprox(refjac), true);

  std::cout << "--> Constructor 5 test ended with success." << std::endl;
}
