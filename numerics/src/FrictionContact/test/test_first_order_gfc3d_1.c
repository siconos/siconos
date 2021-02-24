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
#include "Friction_cst.h"                // for SICONOS_GLOBAL_FRICTION_3D_NSGS
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int n_solvers = 11;
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  for(int d =0; d <n_data; d++)
  {
    // GFC3D, NSGS, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // GFC3D, NSGS, default values.
    // projection on cone for the internal solver, with default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // GFC3D, VI_EG, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_VI_EG);
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 40000;
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-1;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // GFC3D, VI_FPP, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_VI_FPP);
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 40000;
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-1;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ACLMFP, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ACLMFP);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ADMM
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ADMM
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;
    current++;
  }
  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ADMM
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO;
    current++;
  }
  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ADMM
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO;
    collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =  SICONOS_FRICTION_3D_RESCALING_SCALAR;
    current++;
  }
  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ADMM
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO;
    collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =  SICONOS_FRICTION_3D_RESCALING_BALANCING_M;
    current++;
  }
  for(int d =0; d <n_data; d++)
  {
    // GFC3D, ADMM
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] = SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;
    collection[current].options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO;
    collection[current].options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =  SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT;
    current++;
  }

  return collection;


}
