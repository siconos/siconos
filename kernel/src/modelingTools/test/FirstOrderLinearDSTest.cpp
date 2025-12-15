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
#include "FirstOrderLinearDSTest.hpp"

#include "FirstOrderLinearDS.hpp"
//  #include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearDSTest);

void FirstOrderLinearDSTest::setUp() {
  x0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  b0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*b0)(0) = 4;
  (*b0)(1) = 5;
  (*b0)(2) = 6;

  A0 = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::io::readDenseMatrix("matA0.dat"));
}
void FirstOrderLinearDSTest::tearDown() {}

// constructor from initial state only, alias mode
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS1_alias() {
  std::cout << "--> Test: constructor 1 (alias)." << std::endl;
  auto ds =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0, siconos::algebra::alias_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->x0() == *x0, true);
  (*x0) *= 2.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->x0() == *x0, true);
  ds->resetToInitialState();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->LU_M() == nullptr, true);

  // Test compute functions that should do nothing but that should not fail.
  double time = 1.5;
  ds->computefVector(*ds->x(), time);
  ds->computeJacobianfOver_x(*ds->x(), time);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeMMatrix(time);

  // Rhs stuff
  ds->initRhs(time);  // Call computeRhs and computeJacobianRhsOver_x
  siconos::algebra::SiconosMatrix zero_mat = siconos::algebra::SiconosMatrix::Zero(3, 3);
  siconos::algebra::SiconosVector zero_vec = siconos::algebra::SiconosVector::Zero(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", *ds->rhs() == zero_vec, true);

  // Linear DS with time-invariant coeff.
  ds->setConstantA(*A0, siconos::algebra::alias_t);
  ds->setConstantbVector(*b0, siconos::algebra::alias_t);
  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  ds->setConstantMMatrix(Mref, siconos::algebra::alias_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->bVector() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->MMatrix() == Mref, true);

  (*A0)(0, 0) += 3;
  Mref *= 2.;
  (*b0)(0) += 1.2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->bVector() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->MMatrix() == Mref, true);

  // Call computes that should not modify coeff but should not fail
  time = 12;
  ds->computefVector(*ds->x(), time);
  ds->computeJacobianfOver_x(*ds->x(), time);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeMMatrix(time);
  ds->updatePlugins(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->bVector() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->MMatrix() == Mref, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", ds->isTimeInvariant(), true);

  // Rhs stuff
  ds->initRhs(time);  // Call computeRhs and computeJacobianRhsOver_x
  siconos::algebra::SiconosDenseLUMatrix LUM{Mref};
  auto tmp = *A0 * *x0 + *b0;

  auto ref_rhs = LUM.solve(tmp).eval();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test constr 1 (alias) : ", (*ds->rhs() - ref_rhs).norm() < 1e-15, true);

  const auto& b00 = ds->jacobianRhsOver_x();
  auto ref_block = LUM.solve(*A0).eval();
  siconos::algebra::SiconosVector jacoref;
  auto x_size = ds->dimension();
  jacoref.noalias() = Eigen::Map<const Eigen::VectorXd>(ref_block.data(), ref_block.size());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", jacoref.isApprox(b00, 1e-9), true);

  std::cout << "✅ Constructor 1 (alias mode) test ended with success." << std::endl;
}

void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS1_copy() {
  std::cout << "--> Test: constructor 1 (copy)." << std::endl;
  auto ds =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0, siconos::algebra::copy_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->x0() == *x0, true);
  (*x0) *= 2.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->x0() == *x0, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->LU_M() == nullptr, true);

  // Test compute functions that should do nothing but that should not fail.
  double time = 1.5;
  ds->computefVector(*ds->x(), time);
  ds->computeJacobianfOver_x(*ds->x(), time);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeMMatrix(time);

  // Rhs stuff
  ds->initRhs(time);  // Call computeRhs and computeJacobianRhsOver_x
  siconos::algebra::SiconosMatrix zero_mat = siconos::algebra::SiconosMatrix::Zero(3, 3);
  siconos::algebra::SiconosVector zero_vec = siconos::algebra::SiconosVector::Zero(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", *ds->rhs() == zero_vec, true);

  // Linear DS with time-invariant coeff.
  ds->setConstantA(*A0, siconos::algebra::copy_t);
  ds->setConstantbVector(*b0, siconos::algebra::copy_t);
  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  ds->setConstantMMatrix(Mref, siconos::algebra::copy_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->bVector() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->MMatrix() == Mref, true);

  (*A0)(0, 0) += 3;
  Mref *= 2.;
  (*b0)(0) += 1.2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->A() == *A0, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->bVector() == *b0, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (copy) : ", ds->MMatrix() == Mref, false);
  std::cout << "✅ Constructor 1 (copy mode) test ended with success." << std::endl;
}

void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS2_alias() {
  std::cout << "--> Test: constructor 2 (alias)." << std::endl;
  auto ds = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0, *A0, *b0,
                                                                    siconos::algebra::alias_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->x0() == *x0, true);
  (*x0) *= 2.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->x0() == *x0, true);
  ds->resetToInitialState();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->LU_M() == nullptr, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->bVector() == *b0, true);
  double time = 12;
  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  ds->setConstantMMatrix(Mref, siconos::algebra::alias_t);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->MMatrix() == Mref, true);

  (*A0)(0, 0) += 3;
  Mref *= 2.;
  (*b0)(0) += 1.2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->bVector() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->MMatrix() == Mref, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", ds->isTimeInvariant(), true);

  // Rhs stuff
  ds->initRhs(time);  // Call computeRhs and computeJacobianRhsOver_x

  // Test compute functions that should do nothing but that should not fail.
  ds->computefVector(*ds->x(), time);
  ds->computeJacobianfOver_x(*ds->x(), time);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeMMatrix(time);
  ds->updatePlugins(time);

  auto LUM = Mref.lu();
  auto tmp = *A0 * *x0 + *b0;
  auto ref_rhs = LUM.solve(tmp).eval();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test const 2 (alias) : ", (*ds->rhs() - ref_rhs).norm() < 1e-15, true);

  const auto& b00 = ds->jacobianRhsOver_x();
  auto ref_block = LUM.solve(*A0).eval();
  siconos::algebra::SiconosVector jacoref;
  auto x_size = ds->dimension();
  jacoref.noalias() = Eigen::Map<const Eigen::VectorXd>(ref_block.data(), ref_block.size());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (alias) : ", jacoref.isApprox(b00, 1e-9), true);

  std::cout << "✅ Constructor 2 (alias mode) test ended with success." << std::endl;
}

void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS2_copy() {
  std::cout << "--> Test: constructor 2 (copy)." << std::endl;
  auto ds = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0, *A0, *b0,
                                                                    siconos::algebra::copy_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->x0() == *x0, true);
  (*x0) *= 2.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->x0() == *x0, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->LU_M() == nullptr, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->bVector() == *b0, true);
  double time = 12;
  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  ds->setConstantMMatrix(Mref, siconos::algebra::copy_t);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->MMatrix() == Mref, true);

  (*A0)(0, 0) += 3;
  Mref *= 2.;
  (*b0)(0) += 1.2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->A() == *A0, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->bVector() == *b0, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test const 2 (copy) : ", ds->MMatrix() == Mref, false);
  std::cout << "✅ Constructor 2 (copy mode) test ended with success." << std::endl;
}

// DS with time-dependant coeff. set with user-defined functions
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS_plugins() {
  std::cout << "--> Test: constructor 1." << std::endl;

  auto ds =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0, siconos::algebra::copy_t);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test with user-defined op. : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test with user-defined op. : ", ds->x0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test with user-defined op. : ", ds->LU_M() == nullptr, true);

  ds->setComputeAFunction([](double time, Eigen::Ref<siconos::algebra::MapType> result) {
    result.setConstant(1.);
    result *= time;
  });

  ds->setComputebVectorFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        result.setConstant(1.);
        result *= 2. * time;
      });

  ds->setComputeMMatrixFunction([](double time, Eigen::Ref<siconos::algebra::MapType> result) {
    result.setZero();
    result(0, 0) = time;
    result(1, 1) = 2. * time;
    result(2, 2) = 3. * time;
  });

  auto time = 2.;

  // All in one call ...
  ds->updatePlugins(time);

  // or one at a time
  ds->computeA(time);
  ds->computeb(time);
  ds->computeMMatrix(time);

  siconos::algebra::SiconosVector ref_vec{3};
  ref_vec << 1, 1, 1;
  siconos::algebra::SiconosMatrix ref_mat{3, 3};
  ref_mat.setConstant(1.);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test with user-defined op. : ", ds->bVector() == 2 * time * ref_vec, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test with user-defined op. : ", ds->A() == time * ref_mat,
                               true);

  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref.setZero();
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test with user-defined op. : ", ds->MMatrix() == Mref, true);

  // Rhs ...
  ds->initRhs(time);
  auto LUM = Mref.lu();
  auto ref_rhs = LUM.solve(ds->A() * *x0 + ds->bVector()).eval();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test with user-defined op. : ", (*ds->rhs() - ref_rhs).norm() < 1e-15, true);

  const auto& b00 = ds->jacobianRhsOver_x();
  auto ref_block = LUM.solve(ds->A()).eval();
  siconos::algebra::SiconosVector jacoref;
  auto x_size = ds->dimension();
  jacoref.noalias() = Eigen::Map<const Eigen::VectorXd>(ref_block.data(), ref_block.size());

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test with user-defined op. : ", jacoref.isApprox(b00, 1e-12),
                               true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test with user-defined op.  : ", ds->isTimeInvariant(), false);
  std::cout << "✅ user-defined op test ended with success." << std::endl;
}
