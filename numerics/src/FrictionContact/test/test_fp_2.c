/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "ConvexQP_cst.h"                // for SICONOS_CONVEXQP_PGOC_LINESE...
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_TFP, SIC...
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  *number_of_tests = 12; //n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  int d;

  // ========== FC3D_Example1_SBM.dat ========
  d = 0;
  // Tresca FP
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  current++;

  // Panagiotopoulos FP
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_PFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  current++;

  // Tresca FP, internal = ConvexQP, PG cylinder. Default values for params
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  solver_options_update_internal(collection[current].options, 0, SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER);
  current++;

  // Tresca FP, internal = ConvexQP, PG cylinder.
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);

  solver_options_update_internal(collection[current].options, 0, SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER);

  collection[current].options->internalSolvers[0]->dparam[SICONOS_CONVEXQP_PGOC_RHO] = -1.;
  collection[current].options->internalSolvers[0]->dparam[SICONOS_CONVEXQP_PGOC_RHOMIN] = -1.e-6;
  collection[current].options->internalSolvers[0]->iparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MAX_ITER] = 20;
  current++;

  // ========== Confeti-ex13-4contact-Fc3D-SBM.dat ========
  d = 2;
  // Tresca FP
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  collection[current].will_fail = 1;
  current++;

  // Panagiotopoulos FP. Default values for params
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_PFP);
  current++;

  // Tresca FP, internal = ConvexQP, PG cylinder.
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  solver_options_update_internal(collection[current].options, 0,
                                 SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER);
  collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  collection[current].options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  current++;

  // ========== Confeti-ex03-Fc3D-SBM.dat ========
  d = 5;

  // Tresca FP
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 2000;
  current++;

  // Panagiotopoulos FP
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_PFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 2000;
  current++;

  // Tresca FP
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 2000;
  collection[current].options->iparam[1] = 1; // ???
  current++;

  // ========== BoxesStack1-i100000-32.hdf5.dat ========
  d = 6;
  // Tresca FP, set dp[3]
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  collection[current].options->dparam[3] = 1e4; // ???
  // expected to fail
  collection[current].will_fail = 1;
  current++;

  // ========== OneObject-i100000-499.hdf5.dat ========
  d = 9;
  // Tresca FP, set dp[3]
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_TFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  collection[current].options->dparam[3] = 1e4; // ???
  collection[current].options->iparam[1] = 1; // ???
  // expected to fail
  collection[current].will_fail = 1;
  current++;

  *number_of_tests = current;
  return collection;

}
