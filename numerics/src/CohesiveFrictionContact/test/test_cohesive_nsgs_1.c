// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for calloc
#include <string.h>  // for strstr

#include "CohesiveFrictionContact_options.h"     // for SICONOS_COHESIVE_FRICTION_3D_NSGS
#include "NonSmoothGaussSeidel_options.h"        // for SICONOS_NSGS_SHUFFLE, ...
#include "SolverOptions.h"                       // for SolverOptions, solver_option...
#include "cohesiveFrictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                          // for TestCase

/* Problems for which the rigid (no compliance) problem cannot be solved to a tight
 * tolerance by NSGS in a reasonable number of sweeps: the active set reached first has no
 * solution (incompatibility ~1e-5 along the self-stress modes of W). */
static int is_hyperstatic(const char* filename) {
  return strstr(filename, "diamond_sphere_sphere_pack") != NULL;
}

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 7;
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = calloc(*number_of_tests, sizeof(TestCase));

  int topsolver = SICONOS_COHESIVE_FRICTION_3D_NSGS;
  int current = 0;

  // 0. NSGS, default values (light error with full final check), tol 1e-4.
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
    current++;
  }

  // 1. NSGS, full error evaluation, tol 1e-4.
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
    collection[current].options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
        SICONOS_NSGS_ERROR_EVALUATION_FULL;
    current++;
  }

  // 2. NSGS, full error evaluation, shuffled contacts, tol 1e-4.
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
    collection[current].options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
        SICONOS_NSGS_ERROR_EVALUATION_FULL;
    collection[current].options->iparam[SICONOS_NSGS_SHUFFLE] = SICONOS_NSGS_SHUFFLE_TRUE;
    current++;
  }

  // 3. NSGS with projected over-relaxation (omega = 1.7), tol 1e-8.
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
        SICONOS_NSGS_ERROR_EVALUATION_FULL;
    collection[current].options->iparam[SICONOS_NSGS_RELAXATION] =
        SICONOS_NSGS_RELAXATION_TRUE;
    collection[current].options->dparam[SICONOS_NSGS_RELAXATION_VALUE] = 1.7;
    if (is_hyperstatic(data_collection[d])) collection[current].will_fail = 1;
    current++;
  }

  // 4. NSGS with Anderson acceleration (m = 5), tol 1e-8.
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
        SICONOS_NSGS_ERROR_EVALUATION_FULL;
    collection[current].options->iparam[SICONOS_COHESIVE_FRICTION_IPARAM_ANDERSON_MEMORY] = 5;
    if (is_hyperstatic(data_collection[d])) collection[current].will_fail = 1;
    current++;
  }

  // 5. NSGS with compliance on all blocks (epsilon = 1e-3) and omega = 1.7, tol 1e-8.
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
    collection[current].options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
        SICONOS_NSGS_ERROR_EVALUATION_FULL;
    collection[current].options->iparam[SICONOS_NSGS_RELAXATION] =
        SICONOS_NSGS_RELAXATION_TRUE;
    collection[current].options->dparam[SICONOS_NSGS_RELAXATION_VALUE] = 1.7;
    collection[current].options->iparam[SICONOS_COHESIVE_FRICTION_IPARAM_COMPLIANCE] =
        SICONOS_COHESIVE_FRICTION_COMPLIANCE_ALL;
    collection[current].options->dparam[SICONOS_COHESIVE_FRICTION_DPARAM_COMPLIANCE] = 1e-3;
    current++;
  }

  // 6. NSGS with compliance on the cohesive points only (epsilon = 1e-2), tol 1e-8.
  for (int d = 0; d < n_data; d++) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
    collection[current].options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
        SICONOS_NSGS_ERROR_EVALUATION_FULL;
    collection[current].options->iparam[SICONOS_COHESIVE_FRICTION_IPARAM_COMPLIANCE] =
        SICONOS_COHESIVE_FRICTION_COMPLIANCE_COHESIVE;
    collection[current].options->dparam[SICONOS_COHESIVE_FRICTION_DPARAM_COMPLIANCE] = 1e-2;
    current++;
  }

  return collection;
}
