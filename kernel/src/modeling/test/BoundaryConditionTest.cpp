/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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
#include "BoundaryConditionTest.hpp"

#include "BoundaryCondition.hpp"
#include "HarmonicBC.hpp"
#include "SiconosVector.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(BoundaryConditionTest);

void BoundaryConditionTest::testBuildBoundaryConditionFixed() {
  std::cout << "--> Test: constructor for fixed BC \n";

  // Case 1: build-in indices tab
  siconos::modeling::BoundaryCondition bc{
      siconos::modeling::BoundaryCondition::Indices{1, 3, 5}};

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - fixed bc - 1: ", bc.size() == 3, true);

  siconos::modeling::BoundaryCondition::Indices ind{1, 3, 5};

  auto bcind = bc.velocityIndices();
  for (siconos::algebra::Index i = 0; i < bc.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - fixed bc - 2: ", ind[i] == bcind[i], true);

  for (auto& vel : bc.prescribedVelocity()) {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - fixed bc - 3: ", vel == 0., true);
  }

  // Case 2: use an existing tab for indices
  siconos::modeling::BoundaryCondition bc2{ind};

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - fixed bc - 4: ", bc2.size() == 3, true);

  auto bcind2 = bc2.velocityIndices();
  for (siconos::algebra::Index i = 0; i < bc2.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - fixed bc - 5: ", ind[i] == bcind2[i], true);

  for (auto& vel : bc2.prescribedVelocity()) {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - fixed bc - 6: ", vel == 0., true);
  }

  std::cout << "✅ Constructor for fixed BC test ended with success." << std::endl;
}

void BoundaryConditionTest::testBuildBoundaryConditionWithValues() {
  std::cout << "--> Test: constructor with prescribed velocity values \n";

  // Case 1: build-in indices tab
  siconos::algebra::SiconosVector velo{3};
  velo << 1, 2, 3;

  siconos::modeling::BoundaryCondition bc{
      siconos::modeling::BoundaryCondition::Indices{1, 3, 5}, velo};

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - bc with prescribed velocity - 1: ", bc.size() == 3,
                               true);

  siconos::modeling::BoundaryCondition::Indices ind{1, 3, 5};

  auto bcind = bc.velocityIndices();
  for (siconos::algebra::Index i = 0; i < bc.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(
        "test - bc with prescribed velocity - 2: ", ind[i] == bcind[i], true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test - bc with prescribed velocity - 3: ", velo == bc.prescribedVelocity(), true);

  std::cout << "✅ Constructor with prescribed velocity values test ended with success.\n";
}

void BoundaryConditionTest::testBuildHarmonicBoundaryConditionScalar() {
  std::cout << "--> Test: constructor for harmonic bc, scalar coeffs \n";

  // Case 1: build-in indices tab
  siconos::algebra::SiconosVector velo{3};
  velo << 1, 2, 3;

  double a = 1;
  double b = 2.;
  double omega = 3.;
  double phi = 4.;

  siconos::modeling::HarmonicBC bc{siconos::modeling::BoundaryCondition::Indices{1, 3, 5}, a,
                                   b, omega, phi};

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - harmonic bc, scalar coeffs - 1: ", bc.size() == 3,
                               true);

  siconos::modeling::BoundaryCondition::Indices ind{1, 3, 5};

  double time = 2.3;

  double val = a + b * cos(time * omega + phi);
  siconos::algebra::SiconosVector real_bc{3};
  real_bc.setConstant(val);

  auto bcind = bc.velocityIndices();
  for (siconos::algebra::Index i = 0; i < bc.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - harmonic bc, scalar coeffs - 2: ", ind[i] == bcind[i],
                                 true);

  bc.computePrescribedVelocity(time);

  for (auto velo : bc.prescribedVelocity())
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - harmonic bc, scalar coeffs - 3: ", velo == val, true);

  std::cout << "✅ Constructor for harmonic bc, scalar coeffs, test ended with success.\n";
}

void BoundaryConditionTest::testBuildHarmonicBoundaryConditionVector() {
  std::cout << "--> Test: constructor for harmonic bc, scalar coeffs \n";

  // Case 1: build-in indices tab
  siconos::algebra::SiconosVector velo{3};
  velo << 1, 2, 3;

  siconos::algebra::Index size = 3;
  siconos::algebra::SiconosVector a{size};
  siconos::algebra::SiconosVector b{size};
  siconos::algebra::SiconosVector omega{size};
  siconos::algebra::SiconosVector phi{size};

  a.setRandom();
  b.setRandom();
  omega.setRandom();
  phi.setRandom();

  siconos::modeling::HarmonicBC bc{siconos::modeling::BoundaryCondition::Indices{1, 3, 5}, a,
                                   b, omega, phi};

  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - harmonic bc, vector coeffs - 1: ", bc.size() == 3,
                               true);

  siconos::modeling::BoundaryCondition::Indices ind{1, 3, 5};

  double time = 2.3;

  siconos::algebra::SiconosVector real_bc{3};
  real_bc = a.array() + b.array() * (omega.array() * time + phi.array()).cos();

  auto bcind = bc.velocityIndices();
  for (siconos::algebra::Index i = 0; i < bc.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - harmonic bc, vector coeffs - 2: ", ind[i] == bcind[i],
                                 true);

  bc.computePrescribedVelocity(time);

  auto vel = bc.prescribedVelocity();
  for (siconos::algebra::Index i = 0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(
        "test - harmonic bc, scalar coeffs - 3: ", vel[i] == real_bc[i], true);

  std::cout << "✅ Constructor for harmonic bc, vector coeffs, test ended with success.\n";
}
