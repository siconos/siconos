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
#include "Friction_cst.h"
#include "gfc3d_Solvers.h"

int main(void)
{
  int info = 0 ;
  char filename[50] = "./data/LMGC_GFC3D_CubeH8.hdf5";
  printf("Test on %s\n", filename);

  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  info = gfc3d_setDefaultSolverOptions(options, SICONOS_GLOBAL_FRICTION_3D_NSGS);
  options->dparam[0] = 1e-08;
  options->iparam[0] = 100000;

  info = gfc3d_test_function_hdf5(filename, options);

  solver_options_delete(options);
  free(options);
  printf("\nEnd of test on %s\n",filename);
  return info;
}
