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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_PROX
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  *number_of_tests = 3; //n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  int d;
  // ========== FC3D_Example1_SBM.dat ========
  d = 0;
  // Prox, default

  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_PROX);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  current++;

  // ========== BoxesStack1-i100000-32.hdf5.dat ========
  d = 6;
  // Prox, set alpha, many iter.
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_PROX);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA] = 1e4;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000000;
  current++;

  // ========== OneObject-i100000-499.hdf5.dat ========
  d = 9;
  // Prox, set alpha, many iter.
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_PROX);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  collection[current].options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA] = 1e4;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000000;
  current++;


  *number_of_tests = current;
  return collection;

}
