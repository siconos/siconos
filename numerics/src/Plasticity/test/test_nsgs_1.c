/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>  // for malloc

#include "../test-utils/plasticity_test_utils.h"  // for build_test_collection
#include "NumericsFwd.h"                          // for SolverOptions
#include "Plasticity_options.h"                   // for PLASTICITY_2D_NSGS, etc
#include "SolverOptions.h"                        // for SolverOptions, solver_option...
#include "test_utils.h"                           // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 4;
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  /* NSGS + default internal solver (NSN_GP_HYBRID) */
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].will_fail = 0;
    collection[current].options = solver_options_create(PLASTICITY_2D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
    current++;
  }

  /* NSGS + Projection on cone */
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].will_fail = 0;
    collection[current].options = solver_options_create(PLASTICITY_2D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
    solver_options_update_internal(collection[current].options, 0,
                                   PLASTICITY_2D_ONECONE_ProjectionOnCone);
    current++;
  }

  /* NSGS + Projection on cone with local iteration */
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].will_fail = 0;
    collection[current].options = solver_options_create(PLASTICITY_2D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    solver_options_update_internal(collection[current].options, 0,
                                   PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 10;
    current++;
  }

  /* NSGS + NSN (pure Newton) */
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].will_fail = 0;
    collection[current].options = solver_options_create(PLASTICITY_2D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
    solver_options_update_internal(collection[current].options, 0, PLASTICITY_2D_ONECONE_NSN);
    current++;
  }

  return collection;
}
