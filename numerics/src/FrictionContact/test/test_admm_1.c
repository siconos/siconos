/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_FRICTION_3D_ADMM
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for SolverOptions, solver_option...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 3;
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  for (int d = 0; d < n_data; d++) {
    // rho strat = constant
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT;
    current++;
  }

  for (int d = 0; d < n_data; d++) {
    // rho strat = residual balancing
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;
    current++;
  }

  for (int d = 0; d < n_data; d++) {
    // rho strat = residual balancing
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;
    collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =
        SICONOS_FRICTION_3D_RESCALING_SCALAR;
    current++;
  }

  return collection;
}
