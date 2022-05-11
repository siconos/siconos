/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_NSN_AC_TEST
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for solver_options_create, Solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  *number_of_tests = 9; //n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int solvers[] = {SICONOS_FRICTION_3D_NSN_AC, SICONOS_FRICTION_3D_NSN_AC_TEST,
                   SICONOS_FRICTION_3D_NSN_FB, SICONOS_FRICTION_3D_NSN_NM
                  };
  int current = 0;

  int d;

  // ===== RockPile_tob1.dat =====
  d= 7;
  for(int s=0; s<4; ++s)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(solvers[s]);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 5e-2;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    current++;
  }

  collection[0].will_fail= 2; // (FC3D_NSN_AC, on ./data/RockPile_tob1.dat)  is unstable
  collection[1].will_fail= 1; // (FC3D_NSN_AC_TEST, on ./data/RockPile_tob1.dat)  is expected to fail.


  // ===== KaplasTower-i1061-4.hdf5.dat =====
  d = 8;
  collection[current].filename = data_collection[d];
  collection[current].options = solver_options_create(SICONOS_FRICTION_3D_NSN_AC_TEST);
  collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
  collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  //#ifdef WITH_UMFPACK
  collection[current].will_fail = 1; // (FC3D_NSN_AC_TEST, on ./data/KaplasTower-i1061-4.hdf5.dat)  is expected to fail.
  //#endif
  current++;

  for(int s=0; s<4; ++s)
  {
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(solvers[s]);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 500;
    current++;
  }
  //#ifdef WITH_UMFPACK
  collection[5].will_fail = 2; //(FC3D_NSN_AC, on ./data/KaplasTower-i1061-4.hdf5.dat)  is unstable
  collection[6].will_fail = 1; //(FC3D_NSN_AC_TEST, on ./data/KaplasTower-i1061-4.hdf5.dat)  is expected to fail.
  //#endif
  *number_of_tests = current;
  return collection;

}
