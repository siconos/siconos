// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_FRICTION_3D_NSN_AC
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  *number_of_tests = 7;  // n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int solvers[] = {SICONOS_FRICTION_3D_NSN_AC, SICONOS_FRICTION_3D_NSN_FB,
                   SICONOS_FRICTION_3D_NSN_NM};
  int current = 0;

  int d;

  // ===== BoxesStack1-i100000-32.hdf5.dat =====
  d = 6;
  for (int s = 0; s < 3; ++s) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(solvers[s]);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1500;
    current++;
  }
  // NSN_FB expected to fail
  collection[1].will_fail = 1;

  // ===== Rover4396.dat =====
  d = 10;

  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_NSN_AC);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  current++;

  for (int s = 0; s < 3; ++s) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(solvers[s]);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    current++;
  }
  *number_of_tests = current;
  return collection;
}
