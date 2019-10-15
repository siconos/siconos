/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_NSGS
#include "SolverOptions.h"               // for SICONOS_DPARAM_TOL, SICONOS_...
#include "frictionContact_test_utils.h"  // for build_friction_test, build_t...
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int n_solvers = 1;
  *number_of_tests = n_data * n_solvers;
  TestCase * tests_list = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  // "External" solver parameters
  // -> same values for all tests.
  // The differences between tests are only for internal solvers and input data.
  int topsolver = SICONOS_FRICTION_3D_NSGS;
  int dpos[] = {1, SICONOS_DPARAM_TOL};  // ipos = [number of values in parameters list, indices]
  double dparam[] = {1e-5};
  int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};  // ipos = [number of values in parameters list, indices]
  int iparam[] = {10000};

  int current = 0;
  for(int d =0; d <n_data; d++)
    {
      // Quartic solver, set tol and max iter.
      double internal_dparam[] = {1e-6};
      int internal_iparam[] = {10};
      build_friction_test(data_collection[d],
                 topsolver, dpos, dparam, ipos, iparam,
                 SICONOS_FRICTION_3D_ONECONTACT_QUARTIC, dpos, internal_dparam, ipos, internal_iparam,
                 &tests_list[current++]);
    }

  return tests_list;

}
