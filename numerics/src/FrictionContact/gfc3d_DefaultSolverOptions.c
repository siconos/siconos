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
#include <string.h>
#include <time.h>
#include <float.h>
#include "gfc3d_Solvers.h"
#include "fc3d_Solvers.h"
#include "NonSmoothDrivers.h"

#include "numerics_verbose.h"

int gfc3d_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  int info = -1;
  switch (solverId)
  {
  case SICONOS_GLOBAL_FRICTION_3D_NSGS:
  {
    info =    fc3d_nsgs_setDefaultSolverOptions(options);
    options->solverId = SICONOS_GLOBAL_FRICTION_3D_NSGS;
    options->iparam[7] = 10;
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR:
  {
    info =    gfc3d_nonsmooth_Newton_AlartCurnier_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_NSGS_WR:
  {
    info =    gfc3d_nsgs_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR:
  {
    info =    gfc3d_nsgs_velocity_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_PROX_WR:
  {
    info =    gfc3d_proximal_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_DSFP_WR:
  {
    info =    gfc3d_DeSaxceFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_TFP_WR:
  {
    info =    gfc3d_TrescaFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC:
  {
    info =    gfc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_VI_EG:
  {
    info = gfc3d_VI_ExtraGradient_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_VI_FPP:
  {
    info = gfc3d_VI_FixedPointProjection_setDefaultSolverOptions(options);
    break;
  }
  default:
  {
    numerics_error("gfc3d_setDefaultSolverOptions", "Unknown Solver");

  }
  }

  return info;
}
