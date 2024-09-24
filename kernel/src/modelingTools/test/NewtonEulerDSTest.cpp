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
#include "NewtonEulerDSTest.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NewtonEulerDSTest);

void NewtonEulerDSTest::setUp() {
  q0 << 1, 2, 3, 0., 1, 0, 0;
  q01 << 1, 2, 3, 1, 0, 0, 0.;
  velocity0 << 4., 5, 6, 7, 8, 9;
  inertia(0, 0) = 1;
  inertia(1, 1) = 2;
  inertia(2, 2) = 3;
}

void NewtonEulerDSTest::tearDown() {}

// constructor from data
void NewtonEulerDSTest::testBuildNewtonEulerDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;

  auto ds = std::make_shared<siconos::modeling::NewtonEulerDS>(q0, velocity0, mass, inertia);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testBuildNewtonEulerDS1A : ", Type::value(*ds) == Type::NewtonEulerDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1B : ", ds->number() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", ds->dimension() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", ds->getqDim() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", ds->scalarMass() == mass, true);

  siconos::algebra::SiconosMatrix massMatrix{6, 6};
  massMatrix(0, 0) = mass;
  massMatrix(1, 1) = mass;
  massMatrix(2, 2) = mass;
  massMatrix.block<3, 3>(3, 3) = inertia;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", *(ds->mass()) == massMatrix,
                               true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildNewtonEulerDS1D : ", ds->computeKineticEnergy() == 595.0, true);

  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}
// constructor from data
void NewtonEulerDSTest::testNewtonEulerDSQuaternion() {
  std::cout << "--> Test: quaternion 1 from position" << std::endl;

  siconos::algebra::SiconosVector axis{3};
  siconos::algebra::SiconosVector axisref(3);
  double angle = 1e24;
  double angleref = 1e24;

  angle = siconos::geometry::axisAngleFromConfiguration(q0, axis);
  std::cout << "q0 angle : " << angle << std::endl;
  std::cout << "q0 axis : " << std::endl;
  std::cout << axis << "\n";
  axisref(0) = 1.0;
  axisref(1) = 0.0;
  axisref(2) = 0.0;
  angleref = M_PI;

  auto R = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);

  siconos::algebra::SiconosMatrix Rref = Eigen::MatrixXd::Identity(3, 3);

  siconos::geometry::computeRotationMatrix(q0, *R);
  R->display();

  Rref(1, 1) = -1.0;
  Rref(2, 2) = -1.0;

  auto v = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto vref = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*v)(0) = 1.0;
  (*v)(1) = 1.0;
  (*v)(2) = 1.0;

  siconos::geometry::quaternionRotate(q0, *v);
  std::cout << "v : " << std::endl;
  v->display();
  (*vref)(0) = 1.0;
  (*vref)(1) = -1.0;
  (*vref)(2) = -1.0;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionA : ", angle == angleref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionB : ", axis == axisref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionC : ", *(R) == Rref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionC : ", *(v) == *(vref), true);

  angle = siconos::geometry::axisAngleFromConfiguration(q01, axis);
  std::cout << "q01 angle : " << angle << std::endl;
  std::cout << "q01 axis : " << axis << "\n";
  axisref(0) = 0.0;
  axisref(1) = 0.0;
  axisref(2) = 0.0;
  angleref = 0.0;

  Rref.eye();
  siconos::geometry::computeRotationMatrix(q01, *R);
  R->display();
  Rref.display();

  siconos::geometry::quaternionRotate(q01, *v);
  std::cout << "v : " << std::endl;
  v->display();
  (*vref)(0) = 1.0;
  (*vref)(1) = -1.0;
  (*vref)(2) = -1.0;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternion : ", angle == angleref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternion : ", axis == axisref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternion : ", *(R) == Rref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternion : ", *(v) == *(vref), true);

  std::cout << " ---------- test with q02 (rotation of pi/2 about the y-axis)" << std::endl;

  auto q02 = std::make_shared<siconos::algebra::SiconosVector>(7);
  q02->zero();
  angle = M_PI_2;
  axis.setZero();
  axis(1) = 1.0;
  std::cout << "q02 angle : " << angle << std::endl;
  std::cout << "q02 axis : " << axis << "\n";
  siconos::geometry::quaternionFromAxisAngle(axis, angle, *q02);
  std::cout << "q02  : " << std::endl;
  q02->display();

  auto q02ref = std::make_shared<siconos::algebra::SiconosVector>(7);
  q02ref->zero();
  (*q02ref)(3) = cos(angle / 2.0);
  (*q02ref)(5) = sin(angle / 2.0);
  q02ref->display();
  siconos::geometry::computeRotationMatrix(*q02, *R);
  R->display();
  Rref.zero();
  Rref(2, 0) = -1.0;
  Rref(1, 1) = 1.0;
  Rref(0, 2) = 1.0;
  Rref.display();

  siconos::geometry::quaternionRotate(*q02, *v);
  std::cout << "v : " << std::endl;
  v->display();
  (*vref)(0) = -1.0;
  (*vref)(1) = -1.0;
  (*vref)(2) = -1.0;
  std::cout << "vref : " << std::endl;
  vref->display();
  auto diff = (*v) - (*vref);
  std::cout << diff.normInf() << std::endl;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionG : ", *(q02) == *(q02ref), true);
  R->display();
  Rref.display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionH : ", R->isApprox(Rref), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternion : ", v->isApprox(*vref), true);

  std::cout << " ---------- test with q03 (rotation of pi/4 about the y-axis)" << std::endl;

  auto q03 = std::make_shared<siconos::algebra::SiconosVector>(7);
  q03->zero();
  angle = M_PI_4;
  axis.setZero();
  axis(1) = 1.0;
  std::cout << "q03 angle : " << angle << std::endl;
  std::cout << "q03 axis : " << axis << "\n";
  siconos::geometry::quaternionFromAxisAngle(axis, angle, *q03);
  std::cout << "q03  : " << std::endl;
  q03->display();

  auto q03ref = std::make_shared<siconos::algebra::SiconosVector>(7);
  q03ref->zero();
  (*q03ref)(3) = cos(angle / 2.0);
  (*q03ref)(5) = sin(angle / 2.0);
  q03ref->display();
  siconos::geometry::computeRotationMatrix(*q03, *R);
  R->display();
  Rref.zero();
  Rref(0, 0) = sqrt(2.0) / 2.0;
  Rref(0, 2) = sqrt(2.0) / 2.0;

  Rref(1, 1) = 1.0;
  Rref(2, 0) = -sqrt(2.0) / 2.0;
  Rref(2, 2) = sqrt(2.0) / 2.0;
  Rref.display();

  siconos::geometry::quaternionRotate(*q03, *v);
  std::cout << "v : " << std::endl;
  v->display();
  (*vref)(0) = -sqrt(2.0);
  (*vref)(1) = -1.0;
  (*vref)(2) = 0.0;
  std::cout << "vref : " << std::endl;
  vref->display();
  auto diff2 = (*v) - (*vref);
  std::cout << diff2.normInf() << std::endl;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionG : ", *(q03) == *(q03ref), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNewtonEulerDSQuaternionH : ", R->isApprox(Rref), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testNewtonEulerDSQuaternion : ",
      (diff2.normInf() <= std::numeric_limits<double>::epsilon() * 10.0), true);

  std::cout << "--> quaternion 2 test ended with success." << std::endl;
}

void NewtonEulerDSTest::testNewtonEulerDSQuaternionMatrix() {
  std::cout << "--> Test: quaternion 2" << std::endl;
  std::cout << " ---------- test with q03 (rotation of pi/4 about the y-axis)" << std::endl;
  auto q03 = std::make_shared<siconos::algebra::SiconosVector>(7);
  q03->zero();
  double angle = M_PI_4;
  siconos::algebra::SiconosVector axis{3};

  axis.setZero();
  axis(1) = 1.0;
  std::cout << "q03 angle : " << angle << std::endl;
  std::cout << "q03 axis : " << axis << "\n";
  siconos::geometry::quaternionFromAxisAngle(axis, angle, *q03);
  std::cout << "q03  : " << std::endl;
  q03->display();

  auto v = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto vref = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*v)(0) = 1.0;
  (*v)(1) = 1.0;
  (*v)(2) = 1.0;
  (*vref)(0) = sqrt(2.0);
  (*vref)(1) = 1.0;
  (*vref)(2) = 0.0;
  std::cout << "v : " << std::endl;
  v->display();
  std::cout << "vref : " << std::endl;
  vref->display();
  siconos::algebra::SiconosVector diff = (*v) - (*vref);
  std::cout << diff.normInf() << std::endl;

  // Old version
  siconos::algebra::SiconosVector aux(3);
  auto matrix = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  siconos::geometry::computeRotationMatrix(*q03, *matrix);  // compute R
  siconos::algebra::prod(*matrix, *v, aux);                 // multiply by R
  *v = aux;
  std::cout << "v : " << std::endl;
  v->display();
  diff = *v - *vref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testNewtonEulerDSQuaternion : ",
      (diff.normInf() <= std::numeric_limits<double>::epsilon() * 10.0), true);

  // double transpose version
  (*v)(0) = 1.0;
  (*v)(1) = 1.0;
  (*v)(2) = 1.0;
  siconos::geometry::computeRotationMatrixTransposed(*q03,
                                                     *matrix);  // Compute R^T for the moment
  siconos::algebra::prod(*v, *matrix, aux);                     // multiply by R^T^T
  *v = aux;
  std::cout << "v : " << std::endl;
  v->display();
  diff = *v - *vref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testNewtonEulerDSQuaternion : ",
      (diff.normInf() <= std::numeric_limits<double>::epsilon() * 10.0), true);

  // New version
  (*v)(0) = 1.0;
  (*v)(1) = 1.0;
  (*v)(2) = 1.0;
  siconos::geometry::quaternionRotate(*q03, *v);
  std::cout << "v : " << std::endl;
  v->display();
  diff = *v - *vref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testNewtonEulerDSQuaternion : ",
      (diff.normInf() <= std::numeric_limits<double>::epsilon() * 10.0), true);

  auto m = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  m->zero();
  (*m)(2, 0) = 1.0;
  (*m)(0, 1) = 1.0;
  (*m)(0, 2) = 1.0;
  (*m)(1, 2) = 1.0;
  (*m)(2, 2) = 1.0;
  auto mref = std::make_shared<siconos::algebra::SiconosMatrix>(3, 3);
  mref->zero();
  (*mref)(0, 0) = sqrt(2.0) / 2.0;
  (*mref)(2, 0) = sqrt(2.0) / 2.0;
  (*mref)(0, 1) = sqrt(2.0) / 2.0;
  (*mref)(2, 1) = -sqrt(2.0) / 2.0;
  (*mref)(0, 2) = sqrt(2.0);
  (*mref)(1, 2) = 1.0;

  siconos::geometry::quaternionRotate(*q03, *m);
  std::cout << "m : " << std::endl;
  m->display();
  std::cout << "mref : " << std::endl;
  mref->display();
  siconos::algebra::SiconosMatrix diff2 = *m - *mref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testNewtonEulerDSQuaternion : ",
      (diff2.normInf() <= std::numeric_limits<double>::epsilon() * 10.0), true);
}

// void NewtonEulerDSTest::testcomputeDS()
// {
//   std::cout << "-->Test: computeDS." <<std::endl;
//   auto  ds = std::make_shared<NewtonEulerDS>(tmpxml2));
//   std::shared_ptr<NewtonEulerDS> copy =  std::static_pointer_cast<NewtonEulerDS>(ds);
//   double time = 1.5;
//   ds->initialize("EventDriven", time);
//   ds->computeRhs(time);
//   std::cout << "-->Test: computeDS." <<std::endl;
//   ds->computeJacobianRhsx(time);
//   std::cout << "-->Test: computeDS." <<std::endl;
//   SiconosMatrix M(3, 3);
//   M(0, 0) = 1;
//   M(1, 1) = 2;
//   M(2, 2) = 3;
//   auto jx = ds->jacobianRhsx();
//   auto vf = ds->rhs();

//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSI : ", *(vf->vector(0)) == *velocity0, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSJ : ", siconos::algebra::prod(M,
//   *(vf->vector(1))) == (copy->getFExt() - copy->getFInt() - copy->getFGyr()) , true);

//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", siconos::algebra::prod(M, *(jx->block(1,
//   0))) == (copy->getJacobianFL(0)) , true); CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL :
//   ", siconos::algebra::prod(M, *(jx->block(1, 1))) == (copy->getJacobianFL(1)) , true);
//   std::cout << "--> computeDS test ended with success." <<std::endl;

// }
