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
#include "SOCLCP_Solvers.h"
#include "NonSmoothDrivers.h"

int soclcp_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  options->iparam = NULL;
  options->dparam = NULL;
  null_SolverOptions(options);

  int info = -1;
  switch(solverId)
  {
  case SICONOS_SOCLCP_NSGS:
  {
    info =    soclcp_nsgs_setDefaultSolverOptions(options);
    break;
  }
  /* case SICONOS_SOCLCP_NSGSV: */
  /* { */
  /*   info =    soclcp_nsgs_velocity_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_*/
  /* { */
  /*   info =    soclcp_proximal_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_TFP: */
  /* { */
  /*   info =    soclcp_TrescaFixedPoint_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_DSFP: */
  /* { */
  /*   info =    soclcp_DeSaxceFixedPoint_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_FPP: */
  /* { */
  /*   info =    soclcp_fixedPointProjection_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_EG: */
  /* { */
  /*   info =    soclcp_ExtraGradient_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  case SICONOS_SOCLCP_VI_FPP:
  {
    info =    soclcp_VI_FixedPointProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_SOCLCP_VI_EG:
  {
    info =    soclcp_VI_ExtraGradient_setDefaultSolverOptions(options);
    break;
  }
  /* case SICONOS_SOCLCP_HP: */
  /* { */
  /*   info =    soclcp_HyperplaneProjection_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_LOCALAC: */
  /* { */
  /*   info =    soclcp_AlartCurnier_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_LOCALFB: */
  /* { */
  /*   info =    soclcp_FischerBurmeister_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_QUARTIC: */
  /* { */
  /*   info =    soclcp_unitary_enumerative_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_QUARTIC_NU: */
  /* { */
  /*   info =    soclcp_unitary_enumerative_setDefaultSolverOptions(options); */
  /*   options->solverId = SICONOS_SOCLCP_QUARTIC_NU; */
  /*   break; */
  /* } */
  default:
  {
    numericsError("soclcp_setDefaultSolverOptions", "Unknown Solver");

  }
  }

  return info;
}
