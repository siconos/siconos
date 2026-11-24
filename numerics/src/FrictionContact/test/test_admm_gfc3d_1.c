// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_GLOBAL_FRICTION_3D_ADMM
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 1;
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, NSGS_WR. */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);
   */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000; */
  /*   current++; */
  /* } */

  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM. */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000000; */
  /*   collection[current].will_fail = 1; // expected to fail */
  /*   current++; */
  /* } */

  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM, set rho strategy */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
   * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].will_fail = 1; // expected to fail */
  /*   current++; */
  /* } */

  for (int d = 0; d < n_data; d++) {
    // GFC3D, ADMM, set rho strategy.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING;
    current++;
  }

  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM, set rho strategy and rescaling */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 30000; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
   * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =
   * SICONOS_FRICTION_3D_RESCALING_SCALAR; */
  /*   current++; */
  /* } */

  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM, set rho strategy and rescaling */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 30000; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
   * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =
   * SICONOS_FRICTION_3D_RESCALING_BALANCING_M; */
  /*   current++; */
  /* } */
  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM, set rho strategy and rescaling */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
   * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =
   * SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT; */
  /*   current++; */
  /* } */
  return collection;
}
