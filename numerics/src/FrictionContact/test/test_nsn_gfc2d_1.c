/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <stdlib.h>  // for malloc

#include "FrictionContact_options.h"     // for SICONOS_GLOBAL_FRICTION_2D_NSN_AC
#include "SolverOptions.h"               // for solver_options_create
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

typedef struct {
  const char* name;
  int formulation;
  int line_search;
  int warm_start;
} Gfc2dNsnOptionPreset;

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int solver = SICONOS_GLOBAL_FRICTION_2D_NSN_AC;

  /* Grid of options used to compare the Alart--Curnier and Jean--Moreau
   * formulations and to find the best combination of line-search / warm-start
   * settings for GFC2D_NSN_AC.
   *
   * We avoid the combination "warm start without line search" because it is
   * prone to cycling on the larger granular-flow problems. */
  Gfc2dNsnOptionPreset presets[] = {
      {"AC_default", SICONOS_FRICTION_2D_NSN_FORMULATION_ALART_CURNIER,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
       SICONOS_FRICTION_2D_NSN_WARM_START_NO},
      {"AC_LSno_WSno", SICONOS_FRICTION_2D_NSN_FORMULATION_ALART_CURNIER,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_NO, SICONOS_FRICTION_2D_NSN_WARM_START_NO},
      {"AC_LSGP_WSY", SICONOS_FRICTION_2D_NSN_FORMULATION_ALART_CURNIER,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
       SICONOS_FRICTION_2D_NSN_WARM_START_YES},
      {"JM_default", SICONOS_FRICTION_2D_NSN_FORMULATION_JEAN_MOREAU,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
       SICONOS_FRICTION_2D_NSN_WARM_START_NO},
      {"JM_LSno_WSno", SICONOS_FRICTION_2D_NSN_FORMULATION_JEAN_MOREAU,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_NO, SICONOS_FRICTION_2D_NSN_WARM_START_NO},
      {"JM_LSGP_WSY", SICONOS_FRICTION_2D_NSN_FORMULATION_JEAN_MOREAU,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
       SICONOS_FRICTION_2D_NSN_WARM_START_YES},
      {"AC_LSGP_WSAUTO", SICONOS_FRICTION_2D_NSN_FORMULATION_ALART_CURNIER,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
       SICONOS_FRICTION_2D_NSN_WARM_START_AUTO},
      {"JM_LSGP_WSAUTO", SICONOS_FRICTION_2D_NSN_FORMULATION_JEAN_MOREAU,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
       SICONOS_FRICTION_2D_NSN_WARM_START_AUTO}};

  int n_presets = (int)(sizeof(presets) / sizeof(presets[0]));

  *number_of_tests = n_data * n_presets;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  for (int d = 0; d < n_data; ++d) {
    for (int p = 0; p < n_presets; ++p) {
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(solver);
      collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
      collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
      collection[current].options->iparam[SICONOS_FRICTION_2D_NSN_FORMULATION] =
          presets[p].formulation;
      collection[current].options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] =
          presets[p].line_search;
      collection[current].options->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_STRATEGY] =
          presets[p].warm_start;
      collection[current].options->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_MAX_ITER] = 50;

      /* The 230-contact granular-flow problem has a very ill-conditioned mass
       * matrix.  On this problem the Newton solver needs both the Goldstein--
       * Price line search and the VI fixed-point warm start to keep the
       * Jacobian factorizable. */
      if (d == 2 && presets[p].warm_start == SICONOS_FRICTION_2D_NSN_WARM_START_NO) {
        collection[current].will_fail = 1;
      } else {
        collection[current].will_fail = 0;
      }

      current++;
    }
  }

  return collection;
}
