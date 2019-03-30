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
#include <string.h>
#include <time.h>
#include <float.h>

#include "rolling_fc3d_Solvers.h"

#include "NonSmoothDrivers.h"
#include "Newton_methods.h"

int rolling_fc3d_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  solver_options_nullify(options);

  int info = -1;
  switch (solverId)
  {
  case SICONOS_ROLLING_FRICTION_3D_NSGS:
  {
    info =    rolling_fc3d_nsgs_setDefaultSolverOptions(options);
    break;
  }
  default:
  {
    solver_options_set(options, solverId);

  }
  }

  return info;
}
