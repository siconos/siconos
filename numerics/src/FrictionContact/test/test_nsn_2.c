/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_FRICTION_3D_NSN_AC_NEW
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  *number_of_tests = 9;  // n_data * n_solvers;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int solvers[] = {SICONOS_FRICTION_3D_NSN_AC, SICONOS_FRICTION_3D_NSN_AC_NEW,
                   SICONOS_FRICTION_3D_NSN_FB, SICONOS_FRICTION_3D_NSN_NM};
  int current = 0;

  int d;

  // ===== RockPile_tob1.dat =====
  d = 7;
  for (int s = 0; s < 4; ++s) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(solvers[s]);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 5e-2;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    current++;
#ifdef WITH_MUMPS
    if (s == 2) {
      //  test is  unstable with mumps (FC3D_NSN_FB, on ./data/RockPile_tob1.dat)
      collection[current - 1].will_fail = 1;
    }

#endif
  }

  collection[0].will_fail = 2;  // (FC3D_NSN_AC, on ./data/RockPile_tob1.dat)  is unstable
  collection[1].will_fail =
      1;  // (FC3D_NSN_AC_TEST, on ./data/RockPile_tob1.dat)  is expected to fail.
  collection[2].will_fail = 1;  // (FC3D_NSN_FB, on ./data/RockPile_tob1.dat)  is expected to
                                // fail with EUCLIDEAN scaling.
  // ===== KaplasTower-i1061-4.hdf5.dat =====
  d = 8;
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_NSN_AC_NEW);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  collection[current].will_fail =
      1;  // (FC3D_NSN_AC_TEST, on ./data/KaplasTower-i1061-4.hdf5.dat)  is expected to fail.
  current++;

  for (int s = 0; s < 4; ++s) {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(solvers[s]);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    current++;
  }
  collection[5].will_fail =
      2;  //(FC3D_NSN_AC, on ./data/KaplasTower-i1061-4.hdf5.dat)  is unstable
  collection[6].will_fail =
      1;  //(FC3D_NSN_AC_TEST, on ./data/KaplasTower-i1061-4.hdf5.dat)  is expected to fail.
  *number_of_tests = current;
  return collection;
}
