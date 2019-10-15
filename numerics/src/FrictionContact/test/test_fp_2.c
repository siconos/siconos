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

  *number_of_tests = 12; //n_data * n_solvers;
  TestCase * tests_list = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));
  
  int current = 0;
  
  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Tresca FP
    int dpos[] = {1, SICONOS_DPARAM_TOL};
    double dparam[] = {1e-16};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {100};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
  }
  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Panagiotopoulos FP
    int dpos[] = {1, SICONOS_DPARAM_TOL};
    double dparam[] = {1e-16};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {100};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
  }

  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Tresca FP, internal = ConvexQP, PG cylinder. Default values for params
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, NULL, NULL, NULL, NULL,
               SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder, NULL, NULL, NULL, NULL,
               &tests_list[current++]);
  }
  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Tresca FP, internal = ConvexQP, PG cylinder.
    //
    int internal_dpos[] = {2, 3, 4};
    double internal_dparam[] = {-1., -1.e-6};
    int internal_ipos[] = {1, 2};
    int internal_iparam[] = {20};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, NULL, NULL, NULL, NULL,
               SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder, internal_dpos, internal_dparam,
               internal_ipos, internal_iparam,
               &tests_list[current++]);
  }
  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Tresca FP
    int dpos[] = {1, SICONOS_DPARAM_TOL};
    double dparam[] = {1e-8};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {1000};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
    // expected to fail
    tests_list[current-1].will_fail = 1;
  }
  
  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Panagiotopoulos FP. Default values for params
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_PFP, NULL, NULL, NULL, NULL,
               -1, NULL, NULL, NULL, NULL,
               &tests_list[current++]);
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Tresca FP, internal = ConvexQP, PG cylinder.
    //
    int dpos[] = {1, SICONOS_DPARAM_TOL};
    double dparam[] = {1.e-4};
    double internal_dparam[] = {1.e-6};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {100};
    int internal_iparam[] = {200};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder, dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
  }

  {
    int d = 5;  // Confeti-ex03-Fc3D-SBM.dat
    // Tresca FP
    int dpos[] = {1, SICONOS_DPARAM_TOL};
    double dparam[] = {1e-8};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {2000};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
  }
  {
    int d = 5;  // Confeti-ex03-Fc3D-SBM.dat
    // Panagiotopoulos FP
    int dpos[] = {1, SICONOS_DPARAM_TOL};
    double dparam[] = {1e-8};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {2000};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
  }
  {
    int d = 5;  // Confeti-ex03-Fc3D-SBM.dat
    // Tresca FP, set ip[1]
    int dpos[] = {1, SICONOS_DPARAM_TOL};
    double dparam[] = {1e-8};
    int ipos[] = {2, SICONOS_IPARAM_MAX_ITER, 1};
    int iparam[] = {2000, 1};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
  }
  {
    int d = 6; // BoxesStack1-i100000-32.hdf5.dat
     // Tresca FP, set dp[3]
    int dpos[] = {2, SICONOS_DPARAM_TOL, 3};
    double dparam[] = {1e-8, 1.e4};
    int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};
    int iparam[] = {2000};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
    // expected to fail
    tests_list[current-1].will_fail = 1;
  }
  {
    int d = 9; // OneObject-i100000-499.hdf5.dat
     // Tresca FP, set ip[1] and dp[3]
    int dpos[] = {2, SICONOS_DPARAM_TOL, 3};
    double dparam[] = {1e-8, 1.e4};
    int ipos[] = {2, SICONOS_IPARAM_MAX_ITER, 1};
    int iparam[] = {10000, 1};
    
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_TFP, dpos, dparam, ipos, iparam,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
    // expected to fail
    tests_list[current-1].will_fail = 1;
  }
  *number_of_tests = current;
  return tests_list;

}
