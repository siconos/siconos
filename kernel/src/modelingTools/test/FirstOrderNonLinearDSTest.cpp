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
#include "FirstOrderNonLinearDSTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderNonLinearDSTest);

void FirstOrderNonLinearDSTest::setUp() {
  xnull = std::make_shared<siconos::algebra::SiconosVector>(3);
  x0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  J0 = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matJ0.dat"));
  M = std::make_shared<siconos::algebra::SiconosMatrix>(
      siconos::algebra::readMatrixFromFile("matM.dat"));
}

void FirstOrderNonLinearDSTest::tearDown() {}

// Test 1 - DS with constant coefficients.
void FirstOrderNonLinearDSTest::testBuildFirstOrderNonLinearDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto ds = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(*x0);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", Type::value(*ds) ==
  // Type::FirstOrderNonLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", ds->x0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", ds->LU_M() == nullptr,
                               true);
  double time = 1.5;

  // Test compute functions that should do nothing but that should not fail.
  ds->computefVector(*ds->x(), time);
  ds->computeJacobianfOver_x(*ds->x(), time);
  ds->computeMMatrix(time);

  // Rhs stuff
  ds->initRhs(time);  // Call computeRhs and computeJacobianRhsOver_x
  siconos::algebra::SiconosMatrix zero_mat = siconos::algebra::SiconosMatrix::Zero(3, 3);
  siconos::algebra::SiconosVector zero_vec = siconos::algebra::SiconosVector::Zero(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", *ds->rhs() == zero_vec,
                               true);

  // Now, let's set constant values for all coeff.
  siconos::algebra::SiconosMatrix ref_mat = siconos::algebra::SiconosMatrix::Random(3, 3);
  siconos::algebra::SiconosVector ref_vec = siconos::algebra::SiconosVector::Random(3);
  ds->setConstantfVector(ref_vec);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", ds->fVector() == ref_vec,
                               true);

  ds->setConstantJacobianfOver_x(ref_mat);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderNonLinearDS1 : ", ds->jacobianfOver_x() == ref_mat, true);

  ds->setConstantMMatrix(ref_mat);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", ds->MMatrix() == ref_mat,
                               true);

  ds->initRhs(time);
  auto LU_M_ = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(ref_mat);
  auto ref_rhs = LU_M_->solve(ref_vec);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", *ds->rhs() == ref_rhs,
                               true);

  auto b00 = ds->jacobianRhsOver_x()->block(0, 0);
  auto ref_block = LU_M_->solve(ref_mat);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", *b00 == ref_block, true);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

// copy constructor
void FirstOrderNonLinearDSTest::testBuildFirstOrderNonLinearDS2() {
  std::cout << "--> Test: constructor 2." << std::endl;
  auto source = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(*x0);
  auto ds = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(*source);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2 : ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2 : ", ds->x0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2 : ", ds->LU_M() == nullptr,
                               true);
  double time = 1.5;
  ds->initRhs(time);
  ds->update(time);
  std::cout << "--> Constructor 2 test ended with success." << std::endl;
}

void FirstOrderNonLinearDSTest::testBuildFirstOrderNonLinearDS3() {
  std::cout << "--> Test: constructor with user-defined functions." << std::endl;
  auto ds = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(*x0);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3 : ", Type::value(*ds) ==
  // Type::FirstOrderNonLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3 : ", ds->dimension() == 3,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3 : ", ds->x0() == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3 : ", ds->LU_M() == nullptr,
                               true);
  double time = 1.5;

  // Set user-defined functions for all operators

  // f(x,t) = t.x
  ds->setComputefVectorFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &x, double time,
         Eigen::Ref<siconos::algebra::MapVectorType> result) { result = time * x; });

  ds->setComputeJacobianfOver_xFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &x, double time,
         Eigen::Ref<siconos::algebra::MapType> result) {
        result << 1, 4, 7, 2, 5, 8, 3, 6, 9;
      });

  ds->setComputeMMatrixFunction([](double time, Eigen::Ref<siconos::algebra::MapType> result) {
    result.setZero();
    result(0, 0) = time;
    result(1, 1) = 2. * time;
    result(2, 2) = 3. * time;
  });

  time = 2.;

  ds->computefVector(*ds->x(), time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderNonLinearDS3 : ", ds->fVector() == time * *x0, true);

  ds->computeJacobianfOver_x(*ds->x(), time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderNonLinearDS3 : ", ds->jacobianfOver_x() == *J0, true);

  ds->computeMMatrix(time);
  siconos::algebra::SiconosMatrix Mref(3, 3);
  Mref.setZero();
  Mref(0, 0) = 1. * time;
  Mref(1, 1) = 2. * time;
  Mref(2, 2) = 3. * time;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3 : ", ds->MMatrix() == Mref,
                               true);

  // Rhs ...
  ds->initRhs(time);     // compute rhs and its jacobian
  auto LUM = Mref.lu();  // Eigen::FullPivLU<siconos::algebra::SiconosMatrix>(Mref);
  auto ref_rhs = LUM.solve(ds->fVector()).eval();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderNonLinearDS3 : ", (*ds->rhs() - ref_rhs).norm() < 1e-15, true);

  auto b00 = ds->jacobianRhsOver_x()->block(0, 0);
  auto ref_block = LUM.solve(ds->jacobianfOver_x()).eval();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderNonLinearDS3 : ", (*b00 - ref_block).norm() < 1e-15, true);

  // Last call to check if everything is ok when we repeat compute... functions
  ds->update(2. * time);

  std::cout << "--> Constructor 3 test ended with success." << std::endl;
}

// init
void FirstOrderNonLinearDSTest::testInitMemory() {
  std::cout << "--> Test: initMemory." << std::endl;
  auto ds1 = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(*xnull);
  ds1->initMemory(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem1 : ", ds1->xMemory().size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem3 : ", ds1->rMemory().size() == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem4 : ", ds1->xMemory().nbVectorsInMemory() == 0,
                               true);
  std::cout << "--> initMemory test ended with success." << std::endl;
}

// swap
void FirstOrderNonLinearDSTest::testSwap() {
  std::cout << "--> Test: swap." << std::endl;
  siconos::algebra::SiconosVector x1{x0->size()};
  x1.setConstant(2.4);
  auto ds = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(x1);
  auto x_ds = ds->x();
  auto r_ds = ds->r();

  *x_ds = x1;
  *r_ds = x1;  // Copy
  ds->initMemory(1);
  ds->swapInMemory();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap1 : ", ds->xMemory().getSiconosVector(0) == x1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap3 : ", ds->rMemory().getSiconosVector(0) == x1, true);
  std::cout << "--> swap test ended with success." << std::endl;
}

// plugins: plugins loading is already in testBuildFirstOrderNonLinearDS2
