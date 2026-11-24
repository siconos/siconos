// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, sqrt
#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_GLOBAL_FRICTION_3D_ADMM
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int n_solvers = 1;
  *number_of_tests = n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  for (int d = 0; d < n_data; d++) {
    // GFC3D,IPM
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_IPM_SNM);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
    current++;
  }
  // for(int d =0; d <n_data; d++)
  // {
  //   // GFC3D,IPM
  //   collection[current].filename = data_collection[d];
  //   collection[current].options =
  //   solver_options_create(SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_PROX);
  //   //   collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  //   // collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  //   SolverOptions * nsn_options = collection[current].options->internalSolvers[0];
  //   nsn_options->dparam[0] = 1e-13;
  //   nsn_options->iparam[0] = 200;
  //   nsn_options->internalSolvers[0]->iparam[0] = 50;
  //   //nsn_options->internalSolvers[0]->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] =
  //   SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO; current++;
  // }

  return collection;
}
