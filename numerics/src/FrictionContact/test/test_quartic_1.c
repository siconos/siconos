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
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_ONECONTA...
#include "frictionContact_test_utils.h"  // for build_friction_test, build_t...
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{

  *number_of_tests = 2; //n_data * n_solvers;
  TestCase * tests_list = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));
  
  int current = 0;
  
  {
    int d = 8; // KaplasTower-i1061-4.hdf5.dat
    // Quartic, default
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_ONECONTACT_QUARTIC, NULL, NULL, NULL, NULL,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
  }
  {
    int d = 9; // OneObject-i100000-499.hdf5.dat
    // Quartic, default
    build_friction_test(data_collection[d],
               SICONOS_FRICTION_3D_ONECONTACT_QUARTIC, NULL, NULL, NULL, NULL,
               -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
  }
  *number_of_tests = current;
  return tests_list;

}
