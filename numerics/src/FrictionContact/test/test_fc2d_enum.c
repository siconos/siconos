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

#include <stdio.h>                       // for NULL
#include <stdlib.h>                      // for malloc
#include "Friction_cst.h"                // for SICONOS_GLOBAL_FRICTION_3D_N...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase
#include "SolverOptions.h"

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int solvers[] = {SICONOS_FRICTION_2D_ENUM};
  
  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));

  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));
  
  int current = 0;
  // tol and maxiter used in tests are the same for all solvers.
  for(int s=0;s<n_solvers;++s)
    {
      for(int d =0; d <n_data; d++)
        {
          // default values for all parameters.
          collection[current].filename = data_collection[d];
          collection[current].options = solver_options_create(solvers[s]);
          collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
          collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
          current++;
        }
    }

  *number_of_tests = current;
  
  return collection;
}
