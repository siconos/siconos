/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <stdlib.h>  // for malloc
#include <string.h>  // for strstr

#include "FrictionContact_options.h"     // for SICONOS_GLOBAL_FRICTION_3D_NSGS
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 3;
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  for (int d = 0; d < n_data; d++) {
    // GFC3D, NSGS, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_2D_NSGS);
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    current++;
  }
  for (int d = 0; d < n_data; d++) {
    // GFC3D, NSGS, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_2D_VI_FPP);
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    current++;
  }
  for (int d = 0; d < n_data; d++) {
    // GFC3D, NSGS, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_2D_VI_EG);
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-2;
    current++;
  }
  /* The VI extra-gradient solver does not reach the requested tolerance on
   * several of the small 2D problems and on two of the granular-flow problems. */
  for (int d = 0; d < n_data; ++d) {
    int idx_eg = 2 * n_data + d;
    if (strstr(data_collection[d], "GFC3D_Example1_2d") ||
        strstr(data_collection[d], "GFC3D_OneContact_2d") ||
        strstr(data_collection[d], "granularflowonwall_00069") ||
        strstr(data_collection[d], "granularflowonwall_00230")) {
      collection[idx_eg].will_fail = 1;
    }
  }

  return collection;
}
