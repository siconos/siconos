/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include <stdio.h>                       // for NULL
#include <stdlib.h>                      // for malloc
#include "Friction_cst.h"
#include "GenericMechanical_cst.h"                // for SICONOS_FRICTION_3D_NSN_HYBR...
#include "SolverOptions.h"               // for SICONOS_DPARAM_TOL, SICONOS_...
#include "genericMechanical_test_utils.h"
#include "test_utils.h"                  // for TestCase
#include "SiconosConfig.h" // for HAS_LAPACK_dgesvd

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
#ifdef HAS_LAPACK_dgesvd
  int n_solvers = 10;
#else
  int n_solvers = 9;
#endif
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));


  // "External" solver parameters
  // -> same values for tol and max iter in all tests.
  // is_reduced parameters set in each test.
  // The differences between tests are only for internal solvers and input data.
  int topsolver = SICONOS_GENERIC_MECHANICAL_NSGS;
  int dpos[] = {1, SICONOS_DPARAM_TOL};  // ipos = [number of values in parameters list, indices]
  double dparam[] = {1e-5};
  int ipos[] = {2, SICONOS_IPARAM_MAX_ITER, SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED};

  int current = 0;

  
  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS};
      // internal = fc3d quartic. 
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_QUARTIC, NULL, NULL, NULL, NULL,
                     &collection[current++]);
      // expected to fail
      collection[5].will_fail = 1;  // GMP5.dat
      collection[6].will_fail = 1;  // GMP6.dat
      
    }

#ifdef HAS_LAPACK_dgesvd
  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIES};
      // internal = fc3d quartic.
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_QUARTIC, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }
#endif

  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES};
      // internal = fc3d quartic.
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_QUARTIC, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }

  for(int d =0; d <n_data-3; d++) // warning : ignore GMP 4 to 6 .dat
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_MLCP_LIKE};
      // internal = fc3d quartic.
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_QUARTIC, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS};
      // internal = fc3d nsn gp
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES};
      // internal = fc3d nsn gp
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }
  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS};
      // internal = fc3d nsn
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_NSN, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES};
      // internal = fc3d nsn
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_NSN, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS};
      // internal = fc3d projection on cone with local iterations
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration, NULL, NULL, NULL, NULL,
                     &collection[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      int iparam[] = {10000, SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES};
      // internal = fc3d projection on cone with local iterations
      build_gmp_test(data_collection[d],
                     topsolver, dpos, dparam, ipos, iparam,
                     SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration, NULL, NULL, NULL, NULL,
                     &collection[current++]);
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
