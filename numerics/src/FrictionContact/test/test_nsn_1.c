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


#include <stdio.h>                       // for NULL
#include <stdlib.h>                      // for malloc
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_NSN_AC_TEST
#include "SolverOptions.h"               // for SICONOS_DPARAM_TOL, SICONOS_...
#include "frictionContact_test_utils.h"  // for build_friction_test, build_test_colle...
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int n_solvers = 6;
  *number_of_tests = n_data * n_solvers;
  TestCase * tests_list = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  for(int d =0; d <n_data; d++)
    {
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-5};
      int ipos[] = {2, SICONOS_IPARAM_MAX_ITER, 1};
      int iparam[] = {10000, 1};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_NSN_AC, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-5};
      int ipos[] = {2, SICONOS_IPARAM_MAX_ITER, 1};
      int iparam[] = {10000, 1};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_NSN_AC_TEST, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
#ifndef WITH_MUMPS
      if(d>=4 && d<9) // Capsules tests work only with mumps
        tests_list[current - 1].will_fail = 1;
#endif
    }

  for(int d =0; d <n_data; d++)
    {
      // change tol and max iter, other params = default.
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-3};
      int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
      int iparam[] = {1000, 1};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_NSN_AC, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      // change tol and max iter, other params = default.
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-3};
      int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
      int iparam[] = {1000, 1};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_NSN_AC_TEST, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
#ifndef WITH_MUMPS
      if(d>=4 && d<9) // Capsules tests work only with mumps
        tests_list[current - 1].will_fail = 1;
#endif
    }

  for(int d =0; d <n_data; d++)
    {
      // change tol and max iter, other params = default.
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-3};
      int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
      int iparam[] = {1000, 1};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_NSN_FB, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      // change tol and max iter, other params = default.
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-3};
      int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
      int iparam[] = {1000, 1};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_NSN_NM, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
    }


  *number_of_tests = current;
  return tests_list;

}
