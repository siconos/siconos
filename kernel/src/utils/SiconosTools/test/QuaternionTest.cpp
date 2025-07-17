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
#include "QuaternionTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(QuaternionTest);

void QuaternionTest::setUp() {}

void QuaternionTest::tearDown() {}

void QuaternionTest::testQuaternion() {
  siconos::algebra::SiconosVector3 axis;
  siconos::algebra::SiconosVector3 axisref;
  double angle = 1e24;
  double angleref = 1e24;

  angle = siconos::geometry::axisAngleFromConfiguration(q0, axis);
  axisref(0) = 1.0;
  axisref(1) = 0.0;
  axisref(2) = 0.0;
  angleref = M_PI;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionA : ", angle == angleref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionB : ", axis == axisref, true);

  siconos::algebra::SiconosMatrix33 R;
  siconos::algebra::SiconosMatrix33 Rref;
  Rref.setZero();
  Rref(0, 0) = 1.;
  Rref(1, 1) = -1.0;
  Rref(2, 2) = -1.0;

  siconos::geometry::computeRotationMatrix(q0, R);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionC : ", R.isApprox(Rref), true);

  siconos::algebra::SiconosVector3 v, vref;
  v.setConstant(1.0);
  vref << 1., -1., -1.;

  siconos::geometry::quaternionRotateVector(q0, v);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionC : ", v == vref, true);

  angle = siconos::geometry::axisAngleFromConfiguration(q01, axis);
  axisref(0) = 0.0;
  axisref(1) = 0.0;
  axisref(2) = 0.0;
  angleref = 0.0;

  Rref = siconos::algebra::SiconosMatrix::Identity(3, 3);
  siconos::geometry::computeRotationMatrix(q01, R);
  siconos::geometry::quaternionRotateVector(q01, v);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternion : ", angle == angleref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternion : ", axis == axisref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternion : ", R == Rref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternion : ", v == vref, true);
  std::cout << "✅ test Quaternion 1 passed.\n";
  angle = M_PI_2;
  axis.setZero();
  axis(1) = 1.0;
  siconos::algebra::SiconosVector7 q02;
  q02.setZero();
  siconos::geometry::quaternionFromAxisAngle(axis, angle, q02);
  siconos::algebra::SiconosVector7 qref;
  qref.setZero();
  qref(3) = cos(angle / 2.0);
  qref(5) = sin(angle / 2.0);
  siconos::geometry::computeRotationMatrix(q02, R);
  Rref.setZero();
  Rref(2, 0) = -1.0;
  Rref(1, 1) = 1.0;
  Rref(0, 2) = 1.0;
  siconos::geometry::quaternionRotateVector(q02, v);
  vref(0) = -1.0;
  vref(1) = -1.0;
  vref(2) = -1.0;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionG : ", q02.isApprox(qref), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionH : ", R.isApprox(Rref), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternion : ", v.isApprox(vref), true);
  std::cout << "✅ test Quaternion (rotation pi/2) passed.\n";

  q02.setZero();
  angle = M_PI_4;
  axis.setZero();
  axis(1) = 1.0;
  siconos::geometry::quaternionFromAxisAngle(axis, angle, q02);
  qref.setZero();
  qref(3) = cos(angle / 2.0);
  qref(5) = sin(angle / 2.0);
  siconos::geometry::computeRotationMatrix(q02, R);
  Rref.setZero();
  Rref(0, 0) = sqrt(2.0) / 2.0;
  Rref(0, 2) = sqrt(2.0) / 2.0;
  Rref(1, 1) = 1.0;
  Rref(2, 0) = -sqrt(2.0) / 2.0;
  Rref(2, 2) = sqrt(2.0) / 2.0;

  siconos::geometry::quaternionRotateVector(q02, v);
  vref(0) = -sqrt(2.0);
  vref(1) = -1.0;
  vref(2) = 0.0;
  auto diff = (v - vref).lpNorm<Eigen::Infinity>();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionG : ", q02.isApprox(qref), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testQuaternionH : ", R.isApprox(Rref), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testQuaternion : ", diff <= std::numeric_limits<double>::epsilon() * 10.0, true);
  std::cout << "✅ test Quaternion (rotation pi/4) passed.\n";
}

void QuaternionTest::testQuaternionMatrix() {
  siconos::algebra::SiconosVector7 q03;
  q03.setZero();
  double angle = M_PI_4;
  siconos::algebra::SiconosVector3 axis;

  axis.setZero();
  axis(1) = 1.0;
  siconos::geometry::quaternionFromAxisAngle(axis, angle, q03);

  siconos::algebra::SiconosVector3 v, vref;
  v.setConstant(1.0);
  vref << sqrt(2.), 1., 0.;

  // Old version
  siconos::algebra::SiconosMatrix33 matrix;
  siconos::geometry::computeRotationMatrix(q03, matrix);  // compute R
  v = matrix * v;                                         // multiply by R
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testQuaternion : ",
      (v - vref).lpNorm<Eigen::Infinity>() <= std::numeric_limits<double>::epsilon() * 10.0,
      true);

  // double transpose version
  v << 1., 1., 1.;
  siconos::geometry::computeRotationMatrixTransposed(q03,
                                                     matrix);  // Compute R^T for the moment

  v = matrix.transpose() * v;  // R^T. v
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testQuaternion : ",
      (v - vref).lpNorm<Eigen::Infinity>() <= std::numeric_limits<double>::epsilon() * 10.0,
      true);

  // New version
  v << 1., 1., 1.;
  siconos::geometry::quaternionRotateVector(q03, v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testQuaternion : ",
      (v - vref).lpNorm<Eigen::Infinity>() <= std::numeric_limits<double>::epsilon() * 10.0,
      true);

  siconos::algebra::SiconosMatrix33 m;
  m.setZero();
  m(2, 0) = 1.0;
  m(0, 1) = 1.0;
  m(0, 2) = 1.0;
  m(1, 2) = 1.0;
  m(2, 2) = 1.0;
  siconos::algebra::SiconosMatrix33 mref;
  mref.setZero();
  mref(0, 0) = sqrt(2.0) / 2.0;
  mref(2, 0) = sqrt(2.0) / 2.0;
  mref(0, 1) = sqrt(2.0) / 2.0;
  mref(2, 1) = -sqrt(2.0) / 2.0;
  mref(0, 2) = sqrt(2.0);
  mref(1, 2) = 1.0;

  siconos::geometry::quaternionRotateMatrix(q03, m);
  auto diff = (m - mref).cwiseAbs().rowwise().sum().maxCoeff();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testQuaternion : ", diff <= std::numeric_limits<double>::epsilon() * 10.0, true);
  std::cout << "✅ test quaternion 2 passed.\n";
}