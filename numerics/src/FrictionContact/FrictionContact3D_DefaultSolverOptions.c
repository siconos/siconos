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
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"

int frictionContact3D_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  null_SolverOptions(options);

  int info = -1;
  switch (solverId)
  {
  case SICONOS_FRICTION_3D_NSGS:
  {
    info =    frictionContact3D_nsgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSGSV:
  {
    info =    frictionContact3D_nsgs_velocity_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_PROX:
  {
    info =    frictionContact3D_proximal_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_TFP:
  {
    info =    frictionContact3D_TrescaFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_ACLMFP:
  {
    info =    frictionContact3D_ACLMFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_SOCLCP:
  {
    info =    frictionContact3D_SOCLCP_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_DSFP:
  {
    info =    frictionContact3D_DeSaxceFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_FPP:
  {
    info =    frictionContact3D_fixedPointProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_EG:
  {
    info =    frictionContact3D_ExtraGradient_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_VI_FPP:
  {
    info =    frictionContact3D_VI_FixedPointProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_VI_EG:
  {
    info =    frictionContact3D_VI_ExtraGradient_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_HP:
  {
    info =    frictionContact3D_HyperplaneProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_LOCALAC:
  {
    info =    frictionContact3D_AlartCurnier_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_LOCALFB:
  {
    info =    frictionContact3D_FischerBurmeister_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_QUARTIC:
  {
    info =    frictionContact3D_unitary_enumerative_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_QUARTIC_NU:
  {
    info =    frictionContact3D_unitary_enumerative_setDefaultSolverOptions(options);
    options->solverId = SICONOS_FRICTION_3D_QUARTIC_NU;
    break;
  }
  default:
  {
    set_SolverOptions(options, solverId);

  }
  }

  return info;
}
