// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_FRICTION_3D_NSGS
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for SolverOptions, solver_option...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 1;
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  // "External" solver parameters
  // -> same values for all tests.
  // The differences between tests are only for internal solvers and input data.
  int topsolver = SICONOS_FRICTION_3D_NSGS;

  int current = 0;
  for (int d = 0; d < n_data; d++) {
    // Quartic solver, set tol and max iter.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-6;
    current++;
  }

  return collection;
}
