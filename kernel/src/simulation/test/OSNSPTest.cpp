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
#include "OSNSPTest.hpp"

#include "FrictionContact.hpp"
#include "NumericsSolversNamespace.h"
#include "SolverOptions.h"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(OSNSPTest);

void OSNSPTest::setUp() {}

void OSNSPTest::tearDown() {}

void OSNSPTest::testOSNSBuild_default()
{
  // Build from solver id
  auto problem = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>();

  auto options = problem->numericsSolverOptions();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test solver options : ",
      options->solverId == SICONOS_FRICTION_3D_NSGS, true);
}

void OSNSPTest::testOSNSBuild_solverid()
{
  // Build from solver id
  auto problem = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(
      3, SICONOS_FRICTION_3D_ADMM);

  auto options = problem->numericsSolverOptions();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test solver options : ",
      options->solverId == SICONOS_FRICTION_3D_ADMM, true);
}

void OSNSPTest::testOSNSBuild_options()
{
  // Build from solver id
  std::shared_ptr<SolverOptions> options{
      solver_options_create(SICONOS_FRICTION_3D_ADMM),
      solver_options_delete};
  auto problem = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(3, options);

  auto options_link = problem->numericsSolverOptions();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "test solver options : ",
      options_link->solverId == SICONOS_FRICTION_3D_ADMM, true);
}
