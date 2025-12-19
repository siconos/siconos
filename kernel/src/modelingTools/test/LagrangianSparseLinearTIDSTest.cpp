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
#include "LagrangianSparseLinearTIDSTest.hpp"

#include "LagrangianSparseLinearTIDS.hpp"
#include "SiconosAlgebraAddons.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
#include "TypeName.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianSparseLinearTIDSTest);

void LagrangianSparseLinearTIDSTest::setUp() {
  q0 << 1, 2, 3;
  velocity0 << 4, 5, 6;

  mass = std::make_shared<siconos::algebra::SiconosSparseMatrix>(3, 3);
  mass->reserve(Eigen::VectorXi::Constant(3, 1));
  mass->insert(0, 0) = 1.;
  mass->insert(1, 1) = 2;
  mass->insert(2, 2) = 3;
  mass->makeCompressed();
}

void LagrangianSparseLinearTIDSTest::tearDown() {}

static void checkRhs(siconos::modeling::LagrangianSparseLinearTIDS& ds,
                     const siconos::algebra::SiconosSparseMatrix& mass,
                     const siconos::algebra::SiconosSparseMatrix& stiffness,
                     const siconos::algebra::SiconosSparseMatrix& damping,
                     const std::string& msg_head) {
  const auto& pos = ds.q_read();
  const auto& velocity = ds.velocity_read();
  siconos::algebra::SiconosVector x0;
  siconos::algebra::SiconosVector rhs0;
  siconos::algebra::concatenateVectors(x0, pos, velocity);
  siconos::algebra::SiconosVector3 ref0;
  siconos::algebra::SiconosSparseLUMatrix solver;
  solver.compute(mass);
  assert(solver.info() == Eigen::Success);
  siconos::algebra::SiconosVector forces = ds.fext() - stiffness * pos - damping * velocity;
  ref0 = solver.solve(forces);
  siconos::algebra::concatenateVectors(rhs0, velocity, ref0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg_head, ds.x_size() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg_head, ds.x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg_head, *(ds.rhs()) == rhs0, true);
  // Compute a reference reference jacobian, for comparison
  siconos::algebra::SiconosVector jacoref{36};
  jacoref.setZero();
  auto ndof = ds.dimension();
  auto idx = [&](int row, int col) { return col * (2 * ndof) + row; };
  for (Eigen::Index i = 0; i < ndof; ++i) jacoref(idx(i, ndof + i)) = 1.0;
  siconos::algebra::SiconosSparseMatrix lower_left = solver.solve(-1. * stiffness);
  siconos::algebra::SiconosSparseMatrix lower_right = solver.solve(-1. * damping);
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
// Basic constructor from mass (copy) and initial state only
//
void LagrangianSparseLinearTIDSTest::testBuildLagrangianSparseLinearTIDS_basic() {
  auto ds = std::make_shared<siconos::modeling::LagrangianSparseLinearTIDS>(
      q0, velocity0, *mass, siconos::algebra::copy_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildLagrangianSparseLinearTIDS1: ",
      siconos::types::type_value(*ds) == siconos::modeling::Type::LagrangianSparseLinearTIDS,
      true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 1: ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 2: ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 3: ", ds->velocity0() == velocity0,
                               true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 4: ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 5: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 6: ", ds->p(2) == nullptr, true);
  auto Ecref = 0.5 * (*mass * velocity0).dot(velocity0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Basic constructor - 7: ", ds->computeKineticEnergy() == Ecref, true);

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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - Basic constructor - 10: ", *(ds->rhs()) == rhs0, true);

  siconos::algebra::SiconosVector jacoref{36};
  jacoref.setZero();
  for (unsigned int j = 0; j < 3; ++j) {
    jacoref((3 + j) * 6 + j) = 1.0;
  }
  const auto& jaco_rhs = ds->jacobianRhsOver_x();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - Basic constructor -11: ", jacoref.isApprox(jaco_rhs, 1e-12), true);
  std::cout << "✅ Basic constructor test ended with success." << std::endl;
}

// A DS with all operators (mass, fext, stiffness, damping)
// Operators set as constant vectors or matrices, using the "alias" version for sparse matrices
//  (no copy, shared memory, memory not owned by the DS)
void LagrangianSparseLinearTIDSTest::testBuildLagrangianSparseLinearTIDS_alias() {
  Eigen::Map<siconos::algebra::SiconosSparseMatrix> mass_map(
      mass->rows(), mass->cols(), mass->nonZeros(), mass->outerIndexPtr(),
      mass->innerIndexPtr(), mass->valuePtr());
  auto ds = std::make_shared<siconos::modeling::LagrangianSparseLinearTIDS>(
      q0, velocity0, mass_map, siconos::algebra::alias_t);

  auto stiffness = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);
  auto damping = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);
  Eigen::Map<siconos::algebra::SiconosSparseMatrix> stiffness_map(
      stiffness.rows(), stiffness.cols(), stiffness.nonZeros(), stiffness.outerIndexPtr(),
      stiffness.innerIndexPtr(), stiffness.valuePtr());
  Eigen::Map<siconos::algebra::SiconosSparseMatrix> damping_map(
      damping.rows(), damping.cols(), damping.nonZeros(), damping.outerIndexPtr(),
      damping.innerIndexPtr(), damping.valuePtr());

  ds->setStiffnessMatrix(stiffness_map, siconos::algebra::alias_t);
  ds->setDampingMatrix(damping_map, siconos::algebra::alias_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianSparseTIDS - Alias case: ", ds->q0() == q0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->mass().isApprox(*mass), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianSparseTIDS - Alias case: ",
                               ds->stiffnessMatrix().isApprox(stiffness), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianSparseTIDS - Alias case: ",
                               ds->dampingMatrix().isApprox(damping), true);

  auto oldvalue = mass_map.coeff(0, 0);
  // Modify ref mass --> this should update ds.mass
  mass_map.coeffRef(0, 0) = 99.;
  auto Ecref = 0.5 * velocity0.dot(*mass * velocity0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->computeKineticEnergy() == Ecref, true);
  mass_map.coeffRef(0, 0) = oldvalue;

  siconos::algebra::SiconosVector3 fext;
  fext << 1, 2, 3;
  ds->setConstantFext(fext, siconos::algebra::alias_t);

  siconos::algebra::SiconosVector3 ref;
  ref << 1, 2, 3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->fext().isApprox(ref), true);

  siconos::algebra::SiconosVector3 q;
  q << 1, 2, 3;
  siconos::algebra::SiconosVector3 vel;
  vel << 4, 5, 6;

  auto external_forces = [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result << 0, time, 2. * time;
  };

  double time = 1.5;
  ds->setComputeFextFunction(external_forces);
  ds->computeTotalForces(velocity0, q0, time);

  siconos::algebra::SiconosVector3 x01;
  x01 << 0, 1, 2;
  siconos::algebra::SiconosVector3 forces;
  forces = stiffness * q0 + damping * velocity0 - time * x01;

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case : ", ds->fext() == time * x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->totalForces() == -1. * forces, true);

  ds->computeTotalForces(vel, q, time);
  fext << 0, time, 2. * time;
  forces = fext - stiffness * q - damping * vel;

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Alias case: ", ds->totalForces().isApprox(forces), true);

  ds->initRhs(time);
  checkRhs(*ds, *mass, stiffness, damping, "test -  LagrangianSparseTIDS - Alias case:");

  ds->computeRhs(time);
  checkRhs(*ds, *mass, stiffness, damping, "test -  LagrangianSparseTIDS - Alias case:");

  std::cout << "✅ LagrangianSparseTIDS - Alias test ended with success." << std::endl;
}

// A DS with all operators (mass, K, C
// Operators set as constant vectors or matrices, using the "copy" version for sparse matrices
void LagrangianSparseLinearTIDSTest::testBuildLagrangianSparseLinearTIDS_copy() {
  auto ds = std::make_shared<siconos::modeling::LagrangianSparseLinearTIDS>(
      q0, velocity0, *mass, siconos::algebra::copy_t);

  auto stiffness = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);
  auto damping = siconos::algebra::generateRandomSparseMatrix(3, 3, 4);
  ds->setStiffnessMatrix(stiffness, siconos::algebra::copy_t);
  ds->setDampingMatrix(damping, siconos::algebra::copy_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianSparseTIDS - Copy case: ", ds->q0() == q0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->mass().isApprox(*mass), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->p_read(1).isZero(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->p(0) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->p(2) == nullptr, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianSparseTIDS - Copy case: ",
                               ds->stiffnessMatrix().isApprox(stiffness), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - LagrangianSparseTIDS - Copy case: ",
                               ds->dampingMatrix().isApprox(damping), true);

  auto Ecref = 0.5 * velocity0.dot(*mass * velocity0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->computeKineticEnergy() == Ecref, true);

  siconos::algebra::SiconosVector3 fext;
  fext << 1, 2, 3;
  ds->setConstantFext(fext, siconos::algebra::alias_t);

  siconos::algebra::SiconosVector3 ref;
  ref << 1, 2, 3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->fext().isApprox(ref), true);

  siconos::algebra::SiconosVector3 q;
  q << 1, 2, 3;
  siconos::algebra::SiconosVector3 vel;
  vel << 4, 5, 6;

  auto external_forces = [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result << 0, time, 2. * time;
  };

  double time = 1.5;
  ds->setComputeFextFunction(external_forces);
  ds->computeTotalForces(velocity0, q0, time);

  siconos::algebra::SiconosVector3 x01;
  x01 << 0, 1, 2;
  siconos::algebra::SiconosVector3 forces;
  forces = stiffness * q0 + damping * velocity0 - time * x01;

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case : ", ds->fext() == time * x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->totalForces() == -1. * forces, true);

  ds->computeTotalForces(vel, q, time);
  fext << 0, time, 2. * time;
  forces = fext - stiffness * q - damping * vel;

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - LagrangianSparseTIDS - Copy case: ", ds->totalForces().isApprox(forces), true);

  ds->initRhs(time);
  checkRhs(*ds, *mass, stiffness, damping, "test -  LagrangianSparseTIDS - Copy case:");

  ds->computeRhs(time);
  checkRhs(*ds, *mass, stiffness, damping, "test -  LagrangianSparseTIDS - Copy case:");

  std::cout << "✅ LagrangianSparseTIDS - Copy case test ended with success." << std::endl;
}
