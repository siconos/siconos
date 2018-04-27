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

#include "fc2d_Solvers.h"
#include "NonSmoothDrivers.h"
#include "numerics_verbose.h"


int fc2d_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  int info = -1;
  switch (solverId)
  {
  case SICONOS_FRICTION_2D_NSGS:
  {
    info =    fc2d_sparse_nsgs_setDefaultSolverOptions(options);
    /* info =    fc2d_nsgs_setDefaultSolverOptions(options); */
    break;
  }
  case SICONOS_FRICTION_2D_PGS:
  {
    info =    fc2d_nsgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_2D_CPG:
  {
    info =    fc2d_cpg_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_2D_LATIN:
  {
    info =    fc2d_latin_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_2D_LEMKE:
  {
    info =    fc2d_lexicolemke_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_2D_ENUM:
  {
    info =    fc2d_enum_setDefaultSolverOptions(options);
    break;
  }

  default:
  {
    numerics_error("fc2d_setDefaultSolverOptions", "Unknown Solver");

  }
  }

  return info;
}
