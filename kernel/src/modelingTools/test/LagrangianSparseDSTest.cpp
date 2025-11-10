/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#include "LagrangianSparseDSTest.hpp"

#include "SiconosAlgebraAddons.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "TypeName.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianSparseDSTest);

void LagrangianSparseDSTest::setUp() {
  q0 << 1, 2, 3;
  velocity0 << 4, 5, 6;

  mass = std::make_shared<siconos::algebra::SiconosSparseMatrix>(3, 3);
  mass->reserve(Eigen::VectorXi::Constant(3, 1));
  mass->insert(0, 0) = 1.;
  mass->insert(1, 1) = 2;
  mass->insert(2, 2) = 3;
  mass->makeCompressed();
}

void LagrangianSparseDSTest::tearDown() {}

static void checkRhs(siconos::modeling::LagrangianSparseDS& ds,
                     const siconos::algebra::SiconosSparseMatrix& mass,
                     const siconos::algebra::SiconosSparseMatrix& jacq,
                     const siconos::algebra::SiconosSparseMatrix& jacvel,
                     const siconos::algebra::SiconosVector& ref_forces,
                     const std::string& msg_head) {
  const auto& q0 = ds.q0();
  const auto& velocity0 = ds.velocity0();
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, q0, velocity0);
  siconos::algebra::SiconosVector3 ref0;
  siconos::algebra::SiconosSparseLUMatrix solver;
  solver.compute(mass);
  assert(solver.info() == Eigen::Success);
  ref0 = solver.solve(ref_forces);
  siconos::algebra::concatenateVectors(rhs0, velocity0, ref0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg_head, ds.x_size() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg_head, ds.x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg_head, ds.rhs()->isApprox(rhs0), true);
  // Compute a reference reference jacobian, for comparison
  siconos::algebra::SiconosVector jacoref{36};
  jacoref.setZero();
  auto ndof = ds.dimension();
  auto idx = [&](int row, int col) { return col * (2 * ndof) + row; };
  for (Eigen::Index i = 0; i < ndof; ++i) jacoref(idx(i, ndof + i)) = 1.0;
  siconos::algebra::SiconosSparseMatrix lower_left = solver.solve(-1. * jacq);
  siconos::algebra::SiconosSparseMatrix lower_right = solver.solve(-1. * jacvel);
  for (Eigen::Index k = 0; k < lower_left.outerSize(); ++k) {
    for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(lower_left, k); it; ++it) {
      double value = it.value();
      jacoref(idx(ndof + it.row(), it.col())) = value;
    }
  }
  for (Eigen::Index k = 0; k < lower_left.outerSize(); ++k) {
    for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(lower_right, k); it; ++it) {
      double value = it.value();
      jacoref(idx(ndof + it.row(), ndof + it.col())) = value;
    }
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg_head, ds.jacobianRhsOver_x().isApprox(jacoref), true);
}

// constructor from initial state only
// No plugins, no mass, nothing optional
void LagrangianSparseDSTest::testBuildLagrangianSparseDS_basic() {
  auto ds = std::make_shared<siconos::modeling::LagrangianSparseDS>(q0, velocity0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianSparseDS1: ",
      siconos::types::type_value(*ds) == siconos::modeling::Type::LagrangianSparseDS, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 1: ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 2: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 3: ", ds->velocity0() == velocity0,
                               true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 4: ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 5: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 6: ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Basic constructor - 7: ", ds->computeKineticEnergy() == 38.5, true);

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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 8: ", ds->x_size() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 9: ", ds->x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 10: ", ds->rhs()->isApprox(rhs0),
                               true);

  siconos::algebra::SiconosVector jacoref{36};
  jacoref.setZero();
  for (unsigned int j = 0; j < 3; ++j) {
    jacoref((3 + j) * 6 + j) = 1.0;
  }
  const auto& jaco_rhs = ds->jacobianRhsOver_x();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Basic constructor - 11: ", jacoref.isApprox(jaco_rhs, 1e-12), true);

  ds->computeRhs(time + 0.7);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 12: ", ds->rhs()->isApprox(rhs0),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Basic constructor - 13: ", jacoref.isApprox(jaco_rhs, 1e-12), true);

  std::cout << "✅ Basic constructor test ended with success." << std::endl;
}

// A DS with all operators (mass, fext, fint, fgyr)
// Operators set as constant vectors or matrices, using the "alias" version for sparse matrices
//  (no copy, shared memory, memory not owned by the DS)
void LagrangianSparseDSTest::testBuildLagrangianSparseDS_alias() {
  auto ds = std::make_shared<siconos::modeling::LagrangianSparseDS>(q0, velocity0);
  auto mass_map = std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(
      mass->rows(), mass->cols(), mass->nonZeros(), mass->outerIndexPtr(),
      mass->innerIndexPtr(), mass->valuePtr());

  ds->setConstantMassAlias(*mass_map);

  siconos::algebra::SiconosVector3 fext;
  fext << 1, 2, 3;
  ds->setConstantFext(fext);
  auto oldvalue = mass_map->coeff(0, 0);
  // Modify ref mass --> this should update ds.mass
  mass_map->coeffRef(0, 0) = 99.;

  siconos::algebra::SiconosVector3 ref;
  ref << 1, 2, 3;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Alias case: ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Alias case: ", ds->mass().isApprox(*mass), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ", ds->p_read(1).isZero(),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ", ds->p(0) == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ", ds->p(2) == nullptr,
                               true);
  const auto& velds = ds->velocity_read();
  auto Ecref = 0.5 * velds.dot(*mass * velds);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Alias case: ", ds->computeKineticEnergy() == Ecref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ", ds->fext().isApprox(ref),
                               true);
  mass_map->coeffRef(0, 0) = oldvalue;

  ds->setComputeFintFunction([](const siconos::algebra::SiconosVector& velocity,
                                const siconos::algebra::SiconosVector& q, double time,
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

  // set 'random' value for jacobians, whatever fint is, just for tests
  auto jacq = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);
  auto jacvel = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);
  Eigen::Map<siconos::algebra::SiconosSparseMatrix> jacq_map(
      jacq.rows(), jacq.cols(), jacq.nonZeros(), jacq.outerIndexPtr(), jacq.innerIndexPtr(),
      jacq.valuePtr());
  Eigen::Map<siconos::algebra::SiconosSparseMatrix> jacvel_map(
      jacvel.rows(), jacvel.cols(), jacvel.nonZeros(), jacvel.outerIndexPtr(),
      jacvel.innerIndexPtr(), jacvel.valuePtr());

  ds->setConstantJacobianFintOver_qAlias(jacq_map);
  ds->setConstantJacobianFintOver_velocityAlias(jacvel_map);

  siconos::algebra::SiconosVector3 reffint;
  reffint << 9, 12, 15;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Alias case: ", ds->fint().isApprox(reffint), true);

  ds->computeTotalForces(vel, q, time);
  ds->computeJacobianTotalForcesOver_q(vel, q, time);
  ds->computeJacobianTotalForcesOver_velocity(vel, q, time);

  auto refjac = -1. * jacq;
  siconos::algebra::SiconosVector3 refforces;
  refforces = fext - reffint;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Alias case: ", ds->totalForces().isApprox(refforces), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ",
                               ds->jacobianTotalForcesOver_q().isApprox(refjac), true);
  auto refjac2 = -1. * jacvel;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Alias case: ",
                               ds->jacobianTotalForcesOver_velocity().isApprox(refjac2), true);

  ds->initRhs(time);
  checkRhs(*ds, *mass, jacq, jacvel, refforces, "test -  LagrangianDS - Alias case:");

  time += 0.3;
  ds->computeRhs(time);
  reffint = ds->velocity_read() * time + ds->q_read();
  refforces = fext - reffint;
  checkRhs(*ds, *mass, jacq, jacvel, refforces, "test -  LagrangianDS - Alias case:");

  std::cout << "✅ LagrangianDS - Alias test ended with success." << std::endl;
}

// A DS with all operators (mass, fext, fint, fgyr)
// Operators set as constant vectors or matrices, using the "copy" version for sparse matrices
void LagrangianSparseDSTest::testBuildLagrangianSparseDS_copy() {
  auto ds = std::make_shared<siconos::modeling::LagrangianSparseDS>(q0, velocity0);
  ds->setConstantMassCopy(*mass);

  siconos::algebra::SiconosVector3 fext;
  fext << 1, 2, 3;
  ds->setConstantFext(fext);

  siconos::algebra::SiconosVector3 ref;
  ref << 1, 2, 3;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Copy case: ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ", ds->mass().isApprox(*mass),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ", ds->p_read(1).isZero(),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ", ds->p(2) == nullptr, true);
  const auto& velds = ds->velocity_read();
  auto Ecref = 0.5 * velds.dot(*mass * velds);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Copy case: ", ds->computeKineticEnergy() == Ecref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ", ds->fext().isApprox(ref),
                               true);

  ds->setComputeFintFunction([](const siconos::algebra::SiconosVector& velocity,
                                const siconos::algebra::SiconosVector& q, double time,
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

  // set 'random' value for jacobians, whatever fint is, just for tests
  auto jacq = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);
  auto jacvel = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);

  ds->setConstantJacobianFintOver_qCopy(jacq);
  ds->setConstantJacobianFintOver_velocityCopy(jacvel);

  siconos::algebra::SiconosVector3 reffint;
  reffint << 9, 12, 15;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Copy case: ", ds->fint().isApprox(reffint), true);

  ds->computeTotalForces(vel, q, time);
  ds->computeJacobianTotalForcesOver_q(vel, q, time);
  ds->computeJacobianTotalForcesOver_velocity(vel, q, time);

  auto refjac = -1. * jacq;
  siconos::algebra::SiconosVector3 refforces;
  refforces = fext - reffint;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Copy case: ", ds->totalForces().isApprox(refforces), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ",
                               ds->jacobianTotalForcesOver_q().isApprox(refjac), true);
  auto refjac2 = -1. * jacvel;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Copy case: ",
                               ds->jacobianTotalForcesOver_velocity().isApprox(refjac2), true);

  ds->initRhs(time);
  checkRhs(*ds, *mass, jacq, jacvel, refforces, "test -  LagrangianDS - Copy case:");

  time += 0.3;
  ds->computeRhs(time);
  reffint = ds->velocity_read() * time + ds->q_read();
  refforces = fext - reffint;
  checkRhs(*ds, *mass, jacq, jacvel, refforces, "test -  LagrangianDS - Copy case:");

  std::cout << "✅ LagrangianDS - Copy case test ended with success." << std::endl;
}

// A DS with all operators (mass, fext, fint, fgyr and jacobians) as user-defined functions
void LagrangianSparseDSTest::testBuildLagrangianSparseDS_compute() {
  auto ds = std::make_shared<siconos::modeling::LagrangianSparseDS>(q0, velocity0);

  auto mass_func = [](const Eigen::Ref<const siconos::algebra::SiconosVector>& pos,
                      siconos::algebra::SiconosSparseMatrix& result) {
    // Using coeffRef might be very unefficient, depending on the original structure of the
    // matrix, but this is just for tests.

    result.coeffRef(0, 0) = 1 + 0. * pos(0);  // just to check that pos is callable ...
    result.coeffRef(1, 1) = 2.;
    result.coeffRef(2, 2) = 3.;
    result.makeCompressed();  // Required
    // Note FP: using setFromTriplets might probably be more efficient, in some cases.
  };

  ds->setComputeMassFunction(mass_func);

  ds->setComputeFextFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto i = 0;
        for (auto& v : result) v = time;
      });

  ds->setComputeFgyrFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& vel,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
         Eigen::Ref<siconos::algebra::MapVectorType> result) { result = vel + q; });

  // set 'random' value for jacobians, whatever fgyr is, just for tests
  ds->setComputeJacobianFgyrOver_qFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& v,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
         siconos::algebra::SiconosSparseMatrix& result) {
        result.insert(1, 1) = 2.;
        result.makeCompressed();
      });

  ds->setComputeJacobianFgyrOver_velocityFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& v,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
         siconos::algebra::SiconosSparseMatrix& result) {
        result.insert(1, 1) = 3.;
        result.makeCompressed();
      });

  ds->setComputeFintFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& v,
                                const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                                double time,
                                Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result = v + q;
    result *= time;
  });

  ds->setComputeJacobianFintOver_qFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& v,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
         siconos::algebra::SiconosSparseMatrix& result) {
        result.insert(1, 1) = 4.;
        result.makeCompressed();
      });

  ds->setComputeJacobianFintOver_velocityFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& v,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
         siconos::algebra::SiconosSparseMatrix& result) {
        result.insert(1, 1) = 5.;
        result.makeCompressed();
      });

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Compute case: ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Compute case: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Compute case: ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Compute case: ", ds->p_read(1).isZero(),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Compute case: ", ds->p(0) == nullptr,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Compute case: ", ds->p(2) == nullptr,
                               true);

  siconos::algebra::SiconosVector3 q;
  q << 1, 2, 3;
  siconos::algebra::SiconosVector3 vel;
  vel << 4, 5, 6;
  double time = 2.;

  ds->computeMass(q);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Compute case: ", ds->mass().isApprox(*mass), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Compute case: ", ds->computeKineticEnergy() == 87.0, true);

  ds->computeTotalForces(vel, q, time);
  ds->computeJacobianTotalForcesOver_q(vel, q, time);
  ds->computeJacobianTotalForcesOver_velocity(vel, q, time);

  siconos::algebra::SiconosVector3 refforces;
  refforces.setConstant(time);
  refforces -= (time + 1) * (vel + q);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianDS - Compute case: ", ds->totalForces().isApprox(refforces), true);

  siconos::algebra::SiconosSparseMatrix refjac{3, 3};
  refjac.insert(1, 1) = -6;
  refjac.makeCompressed();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Compute case: ",
                               ds->jacobianTotalForcesOver_q().isApprox(refjac), true);
  refjac.coeffRef(1, 1) = -8;
  refjac.makeCompressed();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianDS - Compute case: ",
                               ds->jacobianTotalForcesOver_velocity().isApprox(refjac), true);

  std::cout << "✅ LagrangianDS - Compute case test ended with success." << std::endl;
}
