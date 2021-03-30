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
#include "Friction_cst.h"                // for SICONOS_ROLLING_FRICTION_3D_...
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for SolverOptions, solver_option...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int n_solvers = 5;
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_ROLLING_FRICTION_3D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 50;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-14;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_ROLLING_FRICTION_3D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone);
    current++;
  }
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_ROLLING_FRICTION_3D_ADMM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    current++;
  }
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_ROLLING_FRICTION_3D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-14;
    collection[current].options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 10;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_ROLLING_FRICTION_3D_NSGS);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-14;
    collection[current].options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 10;
    collection[current].options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION] =
    SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE;
    current++;
  }

  
  *number_of_tests = current;
  return collection;

}
