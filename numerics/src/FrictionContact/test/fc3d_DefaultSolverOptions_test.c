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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"


int main(void)
{
  printf("\n Start of test on Default SolverOptions\n");
  int info = 0 ;
  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));



  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_NSGS);
  solver_options_print(options);
  solver_options_delete(options);

  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_NSGSV);
  solver_options_print(options);
  solver_options_delete(options);

  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_PROX);
  solver_options_print(options);
  solver_options_delete(options);

  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_TFP);
  solver_options_print(options);
  solver_options_delete(options);

  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_DSFP);
  solver_options_print(options);
  solver_options_delete(options);

  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_EG);
  solver_options_print(options);
  solver_options_delete(options);

  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_HP);
  solver_options_print(options);
  solver_options_delete(options);

  free(options);

  printf("\n End of test on Default SolverOptions\n");
  return info;
}
