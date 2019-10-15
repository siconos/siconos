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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_NSN_HYBR...
#include "SolverOptions.h"               // for SICONOS_DPARAM_TOL, SICONOS_...
#include "test_utils.h"                  // for TestCase
#include "frictionContact_test_utils.h"  // for build_friction_test, build_t...

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int solvers[] = {SICONOS_GLOBAL_FRICTION_3D_NSGS_WR, SICONOS_GLOBAL_FRICTION_3D_PROX_WR,
                   SICONOS_GLOBAL_FRICTION_3D_DSFP_WR, SICONOS_GLOBAL_FRICTION_3D_TFP_WR,
                   SICONOS_GLOBAL_FRICTION_3D_ADMM_WR};
  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));

  *number_of_tests = n_data * n_solvers;
  TestCase * tests_list = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));
  
  int current = 0;

  // tol and maxiter used in tests are the same for all solvers.
  int dpos[] = {1, SICONOS_DPARAM_TOL};
  double dparam[] = {1e-5};
  int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
  int iparam[] = {10000}; 


  for(int s=0;s<n_solvers;++s)
    {
      for(int d =0; d <n_data; d++)
        {
          // set tol and maxiter, default values for other parameters.
          build_gfc3d_test(data_collection[d],
                           solvers[s], dpos, dparam, ipos, iparam,
                           -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
        }

    }

  tests_list[7].will_fail = 1; // GFC3D_PROX_WR	./data/GFC3D_Example1.dat 
  tests_list[9].will_fail = 1; // GFC3D_PROX_WR	./data/GFC3D_TwoRods1.dat

  return tests_list;

}
