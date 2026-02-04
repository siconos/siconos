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

#include "NewtonEulerDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NewtonEulerDSTest);

void NewtonEulerDSTest::setUp() {
  q0 << 1, 2, 3, 0., 1, 0, 0;
  q01 << 1, 2, 3, 1, 0, 0, 0.;
  twist0 << 4., 5, 6, 7, 8, 9;
  inertia << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  // inertia(0, 0) = 1;
  // inertia(1, 1) = 2;
  // inertia(2, 2) = 3;
}

void NewtonEulerDSTest::tearDown() {}

// constructor from data
void NewtonEulerDSTest::testBuildNewtonEulerDS1_alias() {
  std::cout << "--> Test: constructor1 (alias)" << std::endl;

  siconos::modeling::NewtonEulerDS neds{q0, twist0, mass, inertia, siconos::algebra::alias_t};

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "test constr 1 (alias)A : ", Type::value(*ds) == Type::NewtonEulerDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias) : ", neds.number() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias)  : ", neds.dimension() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias)  : ", neds.getqDim() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test constr 1 (alias)  : ", neds.scalarMass() == mass, true);

  siconos::algebra::SiconosMatrix massMatrix{6, 6};
  massMatrix(0, 0) = mass;
  massMatrix(1, 1) = mass;
  massMatrix(2, 2) = mass;
  massMatrix.block<3, 3>(3, 3) = inertia;

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildNewtonEulerDS5 : ", neds.totalInertiaMatrix() == massMatrix, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildNewtonEulerDS6 : ", neds.totalInertiaMatrix() == massMatrix, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildNewtonEulerDS7 : ", neds.computeKineticEnergy() == 1921.0, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS8 : ", neds.wrench().isZero(), true);

  auto time = 1.;
  siconos::algebra::SiconosVector twist{6};
  twist = twist0;
  siconos::algebra::SiconosVector q{7};
  q = q0;
  siconos::algebra::SiconosVector mgyr{3};
  auto totalInertiaMatrix = neds.totalInertiaMatrix();
  siconos::modeling::newton_euler::computeMgyr(twist, totalInertiaMatrix, mgyr);
  siconos::algebra::SiconosVector ref{3};
  ref << 454, -908, 454;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS8 : ", mgyr.isApprox(ref), true);
  neds.computeWrench(twist, q, time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS90 : ", neds.wrench().tail(3) == -ref,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS90 : ", neds.wrench().head(3).isZero(),
                               true);

  // Now we try all plugins ...

  neds.setComputeFextFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto i = 0;
        for (auto &v : result) v = time * i++;
      });

  neds.setComputeMextFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto i = 0;
        for (auto &v : result) v = time * i++;
      });

  neds.setComputeFintFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &twist,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto i = 0;
        for (auto &v : result) v = time * i++;
      });

  neds.setComputeJacobianFintOver_qFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &twist,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  neds.setComputeJacobianFintOver_twistFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &twist,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  neds.setComputeMintFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &twist,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto i = 0;
        for (auto &v : result) v = time * i++;
      });

  neds.setComputeJacobianMintOver_qFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &twist,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  neds.setComputeJacobianMintOver_twistFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector> &twist,
         const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
         Eigen::Ref<siconos::algebra::MapType> result) {});

  neds.computeWrench(twist, q, time);
  neds.computeJacobianWrenchOver_twist(twist, q, time);
  neds.computeJacobianWrenchOver_q(twist, q, time);

  std::cout << "--> Constructor test ended with success." << std::endl;
}
// constructor from data
