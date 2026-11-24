// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_FRICTION_3D_PROX
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  *number_of_tests = 3;  // n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

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
