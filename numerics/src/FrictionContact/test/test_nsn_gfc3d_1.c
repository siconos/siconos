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

#include "FrictionContact_options.h"     // for SICONOS_GLOBAL_FRICTION_3D_N...
#include "SolverOptions.h"               // for solver_options_create
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

typedef struct {
  const char* name;
  int formulation;
  int line_search;
  int warm_start;
  int regularization;
  double regul;
} Gfc3dNsnOptionPreset;

TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests) {
  int solver = SICONOS_GLOBAL_FRICTION_3D_NSN_AC;

  /* Grid of options used to compare the Alart--Curnier and Jean--Moreau
   * formulations and to find the best combination of line-search / warm-start
   * settings for GFC3D_NSN_AC.
   *
   * We avoid the combination "warm start without line search" because it is
   * prone to cycling on the larger problems. */
  Gfc3dNsnOptionPreset presets[] = {
      {
          "AC_default",
          SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
          SICONOS_FRICTION_3D_NSN_WARM_START_NO,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {
          "AC_LSno_WSno",
          SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_NO,
          SICONOS_FRICTION_3D_NSN_WARM_START_NO,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {
          "AC_LSGP_WSY",
          SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
          SICONOS_FRICTION_3D_NSN_WARM_START_YES,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {
          "JM_default",
          SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
          SICONOS_FRICTION_3D_NSN_WARM_START_NO,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {
          "JM_LSno_WSno",
          SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_NO,
          SICONOS_FRICTION_3D_NSN_WARM_START_NO,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {
          "JM_LSGP_WSY",
          SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
          SICONOS_FRICTION_3D_NSN_WARM_START_YES,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {
          "AC_LSGP_WSAUTO",
          SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
          SICONOS_FRICTION_3D_NSN_WARM_START_AUTO,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {
          "JM_LSGP_WSAUTO",
          SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD,
          SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
          SICONOS_FRICTION_3D_NSN_WARM_START_AUTO,
          SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO,
          0.0,
      },
      {"AC_LSGP_WSNO_REG", SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD,
       SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE,
       SICONOS_FRICTION_3D_NSN_WARM_START_NO, SICONOS_FRICTION_3D_NSN_REGULARIZATION_YES,
       1e-6}};

  int n_presets = (int)(sizeof(presets) / sizeof(presets[0]));

  *number_of_tests = n_data * n_presets;
  TestCase* collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));

  int current = 0;

  for (int d = 0; d < n_data; ++d) {
    for (int p = 0; p < n_presets; ++p) {
      collection[current].filename = data_collection[d];
      collection[current].options = solver_options_create(solver);
      collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
      collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
      collection[current].options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] =
          presets[p].formulation;
      collection[current].options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] =
          presets[p].line_search;
      collection[current].options->iparam[SICONOS_FRICTION_3D_NSN_WARM_START_STRATEGY] =
          presets[p].warm_start;
      collection[current].options->iparam[SICONOS_FRICTION_3D_NSN_WARM_START_MAX_ITER] = 50;
      collection[current].options->iparam[SICONOS_FRICTION_3D_NSN_REGULARIZATION] =
          presets[p].regularization;
      collection[current].options->dparam[SICONOS_FRICTION_3D_NSN_REGUL] = presets[p].regul;

      collection[current].will_fail = 0;

      current++;
    }
  }

  return collection;
}
