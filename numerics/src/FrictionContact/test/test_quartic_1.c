/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_FRICTION_3D_ONECONTA...
#include "SolverOptions.h"               // for solver_options_create
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  *number_of_tests = 2;  // n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  int d;

  d = 8;  // KaplasTower-i1061-4.hdf5.dat
  // Quartic, default
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
  current++;

  d = 9;  // OneObject-i100000-499.hdf5.dat
  // Quartic, default
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
  current++;

  *number_of_tests = current;
  return collection;
}
