/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
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
#include "Friction_cst.h"                // for SICONOS_GLOBAL_FRICTION_3D_N...
#include "frictionContact_test_utils.h"  // for build_gfc3d_test, build_test...
#include "test_utils.h"                  // for TestCase
#include "SiconosConfig.h" // for WITH_MUMPS // IWYU pragma: keep

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  int solvers[] = {SICONOS_GLOBAL_FRICTION_3D_NSN_AC, SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR};
  
  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));

  *number_of_tests = n_data * n_solvers;
  TestCase * tests_list = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));
  
  int current = 0;

  for(int s=0;s<n_solvers;++s)
    {
      for(int d =0; d <n_data; d++)
        {
          // default values for all parameters.
          build_gfc3d_test(data_collection[d],
                           solvers[s], NULL, NULL, NULL, NULL,
                           -1, NULL, NULL, NULL, NULL, &tests_list[current++]);
        }
    }

  // expected to fail
  tests_list[3].will_fail = 1; // NSN_AC on data/GFC3D_OneContact.dat
# ifndef WITH_MUMPS 
  tests_list[4].will_fail = 1; // NSN_AC on data/GFC3D_TwoRods1.dat, mumps only
#endif 
  
  return tests_list;
}
