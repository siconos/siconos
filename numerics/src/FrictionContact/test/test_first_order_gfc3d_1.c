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
#include "Friction_cst.h"                // for SICONOS_GLOBAL_FRICTION_3D_NSGS
#include "frictionContact_test_utils.h"  // for build_gfc3d_test, build_test...
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
      build_gfc3d_test(data_collection[d],
                          SICONOS_GLOBAL_FRICTION_3D_NSGS, NULL, NULL, NULL, NULL,
                          -1, NULL, NULL, NULL, NULL, &collection[current++]);
    }

  for(int d =0; d <n_data; d++)
    {
      // GFC3D, NSGS, default values.
      // projection on cone for the internal solver, with default values.
      build_gfc3d_test(data_collection[d],
                          SICONOS_GLOBAL_FRICTION_3D_NSGS, NULL, NULL, NULL, NULL,
                          SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone, NULL,  NULL, NULL, NULL,
                          &collection[current++]);
    }
  
  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, VI_EG, default values.
      build_gfc3d_test(data_collection[d],
                          SICONOS_GLOBAL_FRICTION_3D_VI_EG, NULL, NULL, NULL, NULL,
                          -1, NULL, NULL, NULL, NULL, &collection[current++]);
    }

  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, VI_FPP, default values.
      build_gfc3d_test(data_collection[d],
                          SICONOS_GLOBAL_FRICTION_3D_VI_FPP, NULL, NULL, NULL, NULL,
                          -1, NULL, NULL, NULL, NULL, &collection[current++]);
      // Expected to fail
      collection[18].will_fail = 1;  /* GFC3D_VI_FPP	./data/GFC3D_OneContact.dat  */

    }

  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, ACLMFP, default values.
      build_gfc3d_test(data_collection[d],
                          SICONOS_GLOBAL_FRICTION_3D_ACLMFP, NULL, NULL, NULL, NULL,
                          -1, NULL, NULL, NULL, NULL, &collection[current++]);
    }

  for ( int d =0; d <n_data; d++)
    {
      // GFC3D, ADMM, default values.
      build_gfc3d_test(data_collection[d],
                          SICONOS_GLOBAL_FRICTION_3D_ADMM, NULL, NULL, NULL, NULL,
                          -1, NULL, NULL, NULL, NULL, &collection[current++]);
    }


  return collection;

}
