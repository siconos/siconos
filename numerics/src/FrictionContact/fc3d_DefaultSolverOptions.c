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
#include "fc3d_Solvers.h"
#include "NonSmoothDrivers.h"

int fc3d_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  null_SolverOptions(options);

  int info = -1;
  switch (solverId)
  {
  case SICONOS_FRICTION_3D_NSGS:
  {
    info =    fc3d_nsgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSGSV:
  {
    info =    fc3d_nsgs_velocity_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_PROX:
  {
    info =    fc3d_proximal_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_TFP:
  {
    info =    fc3d_TrescaFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_ACLMFP:
  {
    info =    fc3d_ACLMFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_SOCLCP:
  {
    info =    fc3d_SOCLCP_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_DSFP:
  {
    info =    fc3d_DeSaxceFixedPoint_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_FPP:
  {
    info =    fc3d_fixedPointProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_EG:
  {
    info =    fc3d_ExtraGradient_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_VI_FPP:
  {
    info =    fc3d_VI_FixedPointProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_VI_EG:
  {
    info =    fc3d_VI_ExtraGradient_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_HP:
  {
    info =    fc3d_HyperplaneProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_AC:
  {
    info =    fc3d_AlartCurnier_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_FB:
  {
    info =    fc3d_FischerBurmeister_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
  {
    info =    fc3d_unitary_enumerative_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
  {
    info =    fc3d_unitary_enumerative_setDefaultSolverOptions(options);
    options->solverId = SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU;
    break;
  }
  default:
  {
    set_SolverOptions(options, solverId);

  }
  }

  return info;
}
