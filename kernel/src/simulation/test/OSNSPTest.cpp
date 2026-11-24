// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "OSNSPTest.hpp"

#include "FrictionContact.hpp"
#include "NumericsSolversNamespace.h"
#include "SolverOptions.h"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(OSNSPTest);

void OSNSPTest::setUp() {}

void OSNSPTest::tearDown() {}

void OSNSPTest::testOSNSBuild_default() {
  // Build from solver id
  auto problem = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>();

  auto options = problem->numericsSolverOptions();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test solver options : ", options->solverId == SICONOS_FRICTION_3D_NSGS, true);
}

void OSNSPTest::testOSNSBuild_solverid() {
  // Build from solver id
  auto problem = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(
      3, SICONOS_FRICTION_3D_ADMM);

  auto options = problem->numericsSolverOptions();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test solver options : ", options->solverId == SICONOS_FRICTION_3D_ADMM, true);
}

void OSNSPTest::testOSNSBuild_options() {
  // Build from solver id
  std::shared_ptr<SolverOptions> options{solver_options_create(SICONOS_FRICTION_3D_ADMM),
                                         solver_options_delete};
  auto problem =
      std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(3, options);

  auto options_link = problem->numericsSolverOptions();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test solver options : ", options_link->solverId == SICONOS_FRICTION_3D_ADMM, true);
}
