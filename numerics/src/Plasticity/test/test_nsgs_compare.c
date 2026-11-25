// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <math.h>    // for fabs
#include <stdlib.h>  // for malloc

#include "../test-utils/plasticity_test_utils.h"  // for build_test_collection
#include "NumericsFwd.h"                          // for SolverOptions
#include "Plasticity_options.h"                   // for PLASTICITY_2D_NSGS, etc
#include "SiconosBlas.h"
#include "SolverOptions.h"          // for SolverOptions, solver_option...
#include "plasticity_2d_solvers.h"  // for plasticity_2d_nsgs, plasticity_2d_nsgs_generic
#include "test_utils.h"             // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 2; /* Legacy NSGS vs Generic NSGS */
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  /* Legacy NSGS with projection on cone */
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].will_fail = 0;
    collection[current].options = solver_options_create(PLASTICITY_2D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    solver_options_update_internal(collection[current].options, 0,
                                   PLASTICITY_2D_ONECONE_ProjectionOnCone);
    current++;
  }

  /* Generic NSGS with projection on cone */
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].will_fail = 0;
    collection[current].options = solver_options_create(PLASTICITY_2D_NSGS_GENERIC);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    solver_options_update_internal(collection[current].options, 0,
                                   PLASTICITY_2D_ONECONE_ProjectionOnCone);
    current++;
  }

  return collection;
}
