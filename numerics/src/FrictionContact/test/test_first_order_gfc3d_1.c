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
#include <stdlib.h>                      // for malloc
#include "Friction_cst.h"                // for SICONOS_GLOBAL_FRICTION_3D_NSGS
#include "SolverOptions.h"               // for solver_options_create, solve...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int n_solvers = 6;
  *number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));
  
  int current = 0;
  for(int d =0; d <n_data; d++)
    {
      // GFC3D, NSGS, default values.
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
      current++;
    }

  for(int d =0; d <n_data; d++)
    {
      // GFC3D, NSGS, default values.
      // projection on cone for the internal solver, with default values.
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
      solver_options_update_internal(collection[current].options, 0,
                                     SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
      current++;
    }
  
  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, VI_EG, default values.
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_VI_EG);
      current++;
    }

  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, VI_FPP, default values.
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_VI_FPP);
      current++;
      // Expected to fail
      collection[18].will_fail = 1;  /* GFC3D_VI_FPP	./data/GFC3D_OneContact.dat  */
    }

  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, ACLMFP, default values.
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ACLMFP);
      current++;
    }

  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, ADMM, default values.
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_ADMM);
      current++;
    }


  return collection;

}
