/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include <stdlib.h>                        // for malloc
#include "Friction_cst.h"                  // for SICONOS_FRICTION_3D_ONECON...
#include "GenericMechanical_cst.h"         // for SICONOS_GENERIC_MECHANICAL...
#include "NumericsFwd.h"                   // for SolverOptions
#include "SiconosConfig.h"                 // for HAS_LAPACK_dgesvd
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "genericMechanical_test_utils.h"  // for build_test_collection
#include "test_utils.h"                    // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
#ifdef HAS_LAPACK_dgesvd
  int n_solvers = 10;
#else
  int n_solvers = 9;
#endif
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));


  int topsolver = SICONOS_GENERIC_MECHANICAL_NSGS;
  int current = 0;

  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d quartic.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);

    current++;
    // expected to fail
    collection[5].will_fail = 1;  // GMP5.dat
    collection[6].will_fail = 1;  // GMP6.dat
  }

#ifdef HAS_LAPACK_dgesvd
  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d quartic.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIES;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    current++;
  }
#endif

  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d quartic.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    current++;
  }

  for(int d =0; d <n_data-3; d++) // warning : ignore GMP 4 to 6 .dat
  {
    // internal = fc3d quartic.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_MLCP_LIKE;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d nsn gp
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS;
    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    collection[current].options->internalSolvers[1]->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] = SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND;
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d nsn gp
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES;
    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    collection[current].options->internalSolvers[1]->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] = SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND;
    current++;
  }


  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d nsn
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_NSN);
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d nsn
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_NSN);
    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d projection on cone with local iterations
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[1]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[1]->dparam[SICONOS_DPARAM_TOL] = 1e-12;

    current++;
  }

  for(int d =0; d <n_data; d++)
  {
    // internal = fc3d projection on cone with local iterations
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    collection[current].options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] = SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES;

    solver_options_update_internal(collection[current].options, 1, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[1]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[1]->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    current++;
  }

  *number_of_tests = current;


  // tests expected to fail

#ifdef HAS_LAPACK_dgesvd
  collection[5].will_fail = 1; // GMP5.dat
  collection[6].will_fail = 1; // GMP6.dat
  collection[19].will_fail = 1;
  collection[31].will_fail = 1;

  collection[45].will_fail = 1;
  collection[58].will_fail = 1;
  collection[59].will_fail = 1;
  collection[65].will_fail = 1;
#else
  collection[5].will_fail = 1;
  collection[6].will_fail = 1;
  collection[12].will_fail = 1;
  collection[24].will_fail = 1;
  collection[38].will_fail = 1;
  collection[51].will_fail = 1;
  collection[52].will_fail = 1;
  collection[58].will_fail = 1;

#endif

  return collection;

}
