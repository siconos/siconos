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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_ADMM_IPA...
#include "SolverOptions.h"               // for SICONOS_DPARAM_TOL, SICONOS_...
#include "frictionContact_test_utils.h"  // for build_friction_test, build_t...
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int n_solvers = 2;
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));



  int current = 0;
  for(int d =0; d <n_data; d++)
    {
      // rho strat = constant
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-5};
      int ipos[] = {2, SICONOS_IPARAM_MAX_ITER, SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY};
      int iparam[] = {10000, SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_ADMM, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &collection[current++]);
    }
  for(int d =0; d <n_data; d++)
    {
      // rho strat = residual balancing
      int dpos[] = {1, SICONOS_DPARAM_TOL}; 
      double dparam[] = {1e-5};
      int ipos[] = {2, SICONOS_IPARAM_MAX_ITER, SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY};
      int iparam[] = {10000, SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING};
      // 
      build_friction_test(data_collection[d],
                 SICONOS_FRICTION_3D_ADMM, dpos, dparam, ipos, iparam,
                 -1, NULL, NULL, NULL, NULL, &collection[current++]);
    }

  // Expected to fail
  collection[24].will_fail = 1;   /* FC3D ADMM	./data/Confeti-ex13-Fc3D-SBM.dat */
  
  return collection;

}
