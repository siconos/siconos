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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_ADMM_RHO
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  *number_of_tests = 4; //n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  {
    int d = 0; // FC3D_Example1_SBM.dat
    // DeSaxce FP, rho  = 2.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_DSFP);

    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
    collection[current].options->dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 2.;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
    // expected to fail
    collection[current].will_fail = 1;
    current++;
  }


  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // DeSaxce FP, rho  = 5e-3, default for others.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_DSFP);

    collection[current].options->dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 5e3;
    // expected to fail
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 5;  // Confeti-ex03-Fc3D-SBM.dat
    // DeSaxce FP, rho  = 1e-2, default for others.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_DSFP);

    collection[current].options->dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 1e2;
    // expected to fail
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 6; // BoxesStack1-i100000-32.hdf5.dat
    // DeSaxce FP, rho  = 8e4.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_DSFP);

    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 8.e4;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
    // expected to fail
    collection[current].will_fail = 1;
    current++;
  }
  *number_of_tests = current;
  return collection;

}
