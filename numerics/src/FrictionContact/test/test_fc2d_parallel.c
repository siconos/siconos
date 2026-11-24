// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int solvers[] = {
      SICONOS_FRICTION_2D_NSGS_GRAPH,
      SICONOS_FRICTION_2D_NSGS_GRAPH_OPTI,
      SICONOS_FRICTION_2D_NSGS_GRAPH_PERMUT,
      SICONOS_FRICTION_2D_NSGS_PERMUT,
      SICONOS_FRICTION_2D_NSGS_GRAPH_PERMUT_NOCOPY,
      SICONOS_FRICTION_2D_NSGS_GRAPH_PERMUT_FIRSTTOUCH,
  };

  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));

  *number_of_tests = 2 * n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;
  // tol and maxiter used in tests are the same for all solvers.

  for (int s = 0; s < n_solvers; ++s) {
    for (int d = 0; d < n_data; d++) {
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(solvers[s]);
      collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
      collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
      current++;
    }

    for (int d = 0; d < n_data; d++) {
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(solvers[s]);
      collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
      collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
      collection[current].options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 20;
      current++;
    }
  }

  *number_of_tests = current;

  return collection;
}
