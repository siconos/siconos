/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_ONECONTA...
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for SolverOptions, solver_option...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int n_solvers = 5;
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = malloc((*number_of_tests) * sizeof(TestCase));


  // "External" solver parameters
  // -> same values for all tests.

  // The differences between tests are only for internal solvers and input data.
  int topsolver = SICONOS_FRICTION_3D_NSGS;
  int current = 0;

  // nsgs + default values for internal solver.
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    current++;
  }

  // nonsmooth newton 'damped', Moreau-Jean formulation. Default for other parameters
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;

    solver_options_update_internal(collection[current].options, 0, SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    collection[current].options->internalSolvers[0]->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD;

    current++;
  }

  // Projection on cone with local iteration, set tol and max iter.
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 10;
    current++;
  }

  // nonsmooth Newton 'damped', set hybrid strategy
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    collection[current].options->internalSolvers[0]->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] = SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP;
    current++;
    // notice that this strategy is the default one.
  }

  // nonsmooth Newton 'damped', set hybrid strategy
  for(int d =0; d <n_data; d++)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    collection[current].options->internalSolvers[0]->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] = SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP;
    current++;
  }

  return collection;

}
