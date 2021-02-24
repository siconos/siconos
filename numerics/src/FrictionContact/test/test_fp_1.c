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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_ACLMFP
#include "NumericsFwd.h"                 // for SolverOptions
#include "SOCLCP_cst.h"                  // for SICONOS_SOCLCP_VI_EG, SICONO...
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  *number_of_tests = 5; //n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  int d;
  // ===== FC3D_Example1_SBM.dat =====
  d = 0;
  // ACLM fixed point + VI EG as internal solver.
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ACLMFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;

  solver_options_update_internal(collection[current].options, 0, SICONOS_SOCLCP_VI_EG);
  // expected to fail
  collection[current].will_fail = 1;
  current++;

  // ===== Confeti-ex13-4contact-Fc3D-SBM.dat =====
  d = 2;
  // ACLM fixed point + VI FPP as internal solver.
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ACLMFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;

  solver_options_update_internal(collection[current].options, 0, SICONOS_SOCLCP_VI_FPP);
  // expected to fail
  collection[current].will_fail = 1;
  current++;

  // ===== Confeti-ex03-Fc3D-SBM.dat =====
  d = 5;
  // ACLM fixed point
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ACLMFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  collection[current].options->iparam[1] = 1;

  // expected to fail
  collection[current].will_fail = 1;
  current++;

  // SOCLCP, default for all values.
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_SOCLCP);
  // expected to fail
  collection[current].will_fail = 1;
  current++;

  // ===== BoxesStack1-i100000-32.hdf5.dat =====
  d = 6;
  // ACLM fixed point
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ACLMFP);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;

  // expected to fail
  collection[current].will_fail = 1;
  current++;

  *number_of_tests = current;
  return collection;

}
