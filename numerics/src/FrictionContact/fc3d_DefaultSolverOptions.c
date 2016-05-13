/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

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
    info =    fc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_FB:
  {
    info =    fc3d_nonsmooth_Newton_FischerBurmeister_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_NM:
  {
    info =    fc3d_nonsmooth_Newton_NaturalMap_setDefaultSolverOptions(options);
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
