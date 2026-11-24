/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_FRICTION_3D_ADMM_RHO
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  *number_of_tests = 4;  // n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  {
    int d = 0;  // FC3D_Example1_SBM.dat
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
    int d = 2;  // Confeti-ex13-4contact-Fc3D-SBM.dat
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
    int d = 6;  // BoxesStack1-i100000-32.hdf5.dat
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
