/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "NumericsOptions.h"
#include "gfc3d_Solvers.h"
#include "fc3d_Solvers.h"
#include "NonSmoothDrivers.h"


int gfc3d_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  int info = -1;
  switch (solverId)
  {
  case SICONOS_FRICTION_3D_GLOBAL_NSGS:
  {
    info =    fc3d_nsgs_setDefaultSolverOptions(options);
    options->solverId = SICONOS_FRICTION_3D_GLOBAL_NSGS;
    break;
  }
  case SICONOS_FRICTION_3D_GLOBAL_NSN_AC_WR:
  {
    info =    gfc3d_globalAlartCurnier_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_GLOBAL_NSGS_WR:
  {
    info =    gfc3d_nsgs_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_GLOBAL_NSGSV_WR:
  {
    info =    gfc3d_nsgs_velocity_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_GLOBAL_PROX_WR:
  {
    info =    gfc3d_proximal_wr_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_GLOBAL_DSFP_WR:
  {
    info =    gfc3d_DeSaxceFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_GLOBAL_TFP_WR:
  {
    info =    gfc3d_TrescaFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_GLOBAL_AC:
  {
    info =    gfc3d_AlartCurnier_setDefaultSolverOptions(options);
    break;
  }
  default:
  {
    numericsError("gfc3d_setDefaultSolverOptions", "Unknown Solver");

  }
  }

  return info;
}
