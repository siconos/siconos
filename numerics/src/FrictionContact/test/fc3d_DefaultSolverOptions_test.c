/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "frictionContact_test_utils.h"
#include "Friction_cst.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"

int main(void)
{
  printf("\n Start of test on Default SolverOptions\n");
  int info = 0 ;
  SolverOptions * options = NULL;
  int solvers[] = {SICONOS_FRICTION_3D_NSGS, SICONOS_FRICTION_3D_NSGSV, SICONOS_FRICTION_3D_PROX,
                   SICONOS_FRICTION_3D_TFP, SICONOS_FRICTION_3D_DSFP, SICONOS_FRICTION_3D_EG, SICONOS_FRICTION_3D_HP
                  };


  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));

  for(int s=0; s<n_solvers; ++s)
  {
    options = solver_options_create(solvers[s]);
    solver_options_print(options);
    solver_options_delete(options);
    options = NULL;
  }
  printf("\n End of test on Default SolverOptions\n");
  return info;
}
