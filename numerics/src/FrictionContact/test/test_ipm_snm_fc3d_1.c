/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <float.h>   // for DBL_EPSILON
#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_GLOBAL_FRICTION_3D_ADMM
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
    collection[current].options = solver_options_create(SICONOS_FRICTION_3D_IPM_SNM);
    current++;
  }
  collection[current - 1].will_fail = 1;

  return collection;
}
