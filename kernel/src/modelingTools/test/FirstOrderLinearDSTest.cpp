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

#include <limits>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
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

// constructor from initial state only
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS0() {
  std::cout << "--> Test: constructor 0." << std::endl;
  auto ds = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->x0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->LU_M() == nullptr, true);

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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", *ds->rhs() == zero_vec,
                               true);

  // Linear DS with time-invariant coeff.
  ds->setConstantA(*A0);
  ds->setConstantbVector(*b0);
  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  ds->setConstantMMatrix(Mref);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->bVector() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->MMatrix() == Mref, true);

  // Call computes that should not modify coeff but should not fail
  time = 12;
  ds->computefVector(*ds->x(), time);
  ds->computeJacobianfOver_x(*ds->x(), time);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeMMatrix(time);
  ds->updatePlugins(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->A() == *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->bVector() == *b0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->MMatrix() == Mref, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->isTimeInvariant(), true);

  // Rhs stuff
  ds->initRhs(time);  // Call computeRhs and computeJacobianRhsOver_x
  auto LUM = Mref.lu();
  auto tmp = *A0 * *x0 + *b0;
  auto ref_rhs = LUM.solve(tmp).eval();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderNonLinearDS1 : ", (*ds->rhs() - ref_rhs).norm() < 1e-15, true);

  const auto& b00 = ds->jacobianRhsOver_x();
  auto ref_block = LUM.solve(*A0).eval();
  siconos::algebra::SiconosVector jacoref;
  auto x_size = ds->dimension();
  jacoref.noalias() = Eigen::Map<const Eigen::VectorXd>(ref_block.data(), ref_block.size());
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderNonLinearDS3 : ", jacoref.isApprox(b00, 1e-9), true);

  std::cout << "--> Constructor 0 test ended with success." << std::endl;
}

// DS with time-dependant coeff. set with user-defined functions
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;

  auto ds = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->x0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->LU_M() == nullptr, true);

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
      "testBuildFirstOrderLinearDS1 : ", ds->bVector() == 2 * time * ref_vec, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->A() == time * ref_mat,
                               true);

  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref.setZero();
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->MMatrix() == Mref, true);

  // Rhs ...
  ds->initRhs(time);
  auto LUM = Mref.lu();
  auto ref_rhs = LUM.solve(ds->A() * *x0 + ds->bVector()).eval();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderLinearDS1 : ", (*ds->rhs() - ref_rhs).norm() < 1e-15, true);

  const auto& b00 = ds->jacobianRhsOver_x();
  auto ref_block = LUM.solve(ds->A()).eval();
  siconos::algebra::SiconosVector jacoref;
  auto x_size = ds->dimension();
  jacoref.noalias() = Eigen::Map<const Eigen::VectorXd>(ref_block.data(), ref_block.size());

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", jacoref.isApprox(b00, 1e-12),
                               true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->isTimeInvariant(),
                               false);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}
