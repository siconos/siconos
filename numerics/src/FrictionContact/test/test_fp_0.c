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

  *number_of_tests = 4; //n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));
  
  int current = 0;
  
  {
    int d = 0; // FC3D_Example1_SBM.dat
    // DeSaxce FP, rho  = 2.
    int dpos[] = {2, SICONOS_DPARAM_TOL, SICONOS_FRICTION_3D_ADMM_RHO};
    double dparam[] = {1e-8, 2.};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {100000};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_DSFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &collection[current++]);
    // expected to fail
    collection[current-1].will_fail = 1;
  }

  
  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // DeSaxce FP, rho  = 5e-3, default for others.
    int dpos[] = {1, SICONOS_FRICTION_3D_ADMM_RHO};
    double dparam[] = {5e3};    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_DSFP, dpos, dparam, NULL, NULL,
               -1, NULL, NULL, NULL, NULL, &collection[current++]);
    // expected to fail
    collection[current-1].will_fail = 1;
  }
  
  {
    int d = 5;  // Confeti-ex03-Fc3D-SBM.dat
    // DeSaxce FP, rho  = 1e-2, default for others.
    int dpos[] = {1, SICONOS_FRICTION_3D_ADMM_RHO};
    double dparam[] = {1e2};    
    // 
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_DSFP, dpos, dparam, NULL, NULL,
               -1, NULL, NULL, NULL, NULL, &collection[current++]);
    // expected to fail
    collection[current-1].will_fail = 1;
  }
  
  {
    int d = 6; // BoxesStack1-i100000-32.hdf5.dat
    // DeSaxce FP, rho  = 8e4.
    int dpos[] = {2, SICONOS_DPARAM_TOL, SICONOS_FRICTION_3D_ADMM_RHO};
    double dparam[] = {1.e-3, 8.e4};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {100000};
    // 
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_DSFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &collection[current++]);
    // expected to fail
    collection[current-1].will_fail = 1;
  }
  *number_of_tests = current;
  return collection;

}
