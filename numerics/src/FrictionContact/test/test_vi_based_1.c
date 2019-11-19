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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_EG, SICO...
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "VI_cst.h"                      // for SICONOS_VI_DPARAM_RHO, SICON...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  *number_of_tests = 8; //n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  int d;
  // ========== FC3D_Example1_SBM.dat ========
  d = 0; 
  // rho = -1
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_EG);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->dparam[SICONOS_VI_DPARAM_RHO] = -1.;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  current++;

  // rho = 1
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_EG);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->dparam[SICONOS_VI_DPARAM_RHO] = 1.;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  current++;

  // ========== Confeti-ex13-4contact-Fc3D-SBM.dat ========
  d = 2; 
  // EG, rho = -3e-3
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_EG);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->dparam[SICONOS_VI_DPARAM_RHO] = -3.e-3;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  // expected to fail
  collection[current].will_fail = 1;
  current++;

  // EG, rho = -1
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_EG);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  collection[current].options->dparam[SICONOS_VI_DPARAM_RHO] = -1;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  current++;


  // HP
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_HP);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  // expected to fail
  collection[current].will_fail = 1;
  current++;
  
  // ========== Capsules-i101-404.dat ========
  d = 6; 
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_VI_FPP);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000000;
  collection[current].options->iparam[SICONOS_VI_IPARAM_ACTIVATE_UPDATE] = 1;
  // expected to fail
  collection[current].will_fail = 1;
  current++;

  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_FPP);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
  // expected to fail
  collection[current].will_fail = 1;
  current++;


  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_VI_EG);    
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
  collection[current].options->iparam[SICONOS_VI_IPARAM_ACTIVATE_UPDATE] = 1;
  // expected to fail
  collection[current].will_fail = 1;
  current++;


  *number_of_tests = current;
  return collection;

}
