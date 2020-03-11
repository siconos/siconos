/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
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
#include "Friction_cst.h"                // for SICONOS_GLOBAL_FRICTION_3D_N...
#include "SolverOptions.h"               // for solver_options_create
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

#include "SiconosConfig.h" // for WITH_MUMPS // IWYU pragma: keep

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int solvers[] = {SICONOS_GLOBAL_FRICTION_3D_NSN_AC, SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR};

  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));

  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  for(int s=0; s<n_solvers; ++s)
  {
    for(int d =0; d <n_data; d++)
    {
      // default values for all parameters.
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(solvers[s]);
      current++;
    }
  }
#ifdef TEST_NSN_COLLECTION_1
  // expected to fail
  collection[2].will_fail = 1; // GFC3D_NSN_AC, on ./data/GFC3D_Example00_badly_scaled.dat
  collection[5].will_fail = 1; // GFC3D_NSN_AC, on ./data/GFC3D_TwoRods1.dat
  collection[3].will_fail = 1; // NSN_AC on data/GFC3D_OneContact.dat

#ifndef WITH_MUMPS
  collection[4].will_fail = 1; // NSN_AC on data/GFC3D_TwoRods1.dat, mumps only
#endif
#endif
  
  return collection;
}
