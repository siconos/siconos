/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include <stdlib.h>                      // for malloc
#include "Friction_cst.h"                // for SICONOS_GLOBAL_FRICTION_3D_ADMM
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  int n_solvers = 1;
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, NSGS_WR. */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS_WR); */
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
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].will_fail = 1; // expected to fail */
  /*   current++; */
  /* } */

  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ADMM, set rho strategy.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING;
    current++;
  }

  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM, set rho strategy and rescaling */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 30000; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_SCALAR; */
  /*   current++; */
  /* } */

  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM, set rho strategy and rescaling */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 30000; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_BALANCING_M; */
  /*   current++; */
  /* } */
  /* for(int d =0; d <n_data; d++) */
  /* { */
  /*   // GFC3D, ADMM, set rho strategy and rescaling */
  /*   collection[current].filename = data_collection[d]; */
  /*   collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5; */
  /*   collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING; */
  /*   collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT; */
  /*   current++; */
  /* } */
  return collection;
}
