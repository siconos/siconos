/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "gfc3d_Solvers.h"

int main(void)
{
  int info = 0 ;
  const char filename[50] = "./data/LMGC_GFC3D_CubeH8.hdf5";
  printf("Test on %s\n", filename);

  TestCase current;
  current.filename = filename;
  current.options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
  current.options->dparam[SICONOS_DPARAM_TOL] = 1e-08;
  current.options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;

  info = globalFrictionContact_test_function(&current);

  solver_options_delete(current.options);
  printf("\nEnd of test on %s\n",filename);
  return info;
}
