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
#include <stdio.h>             // for NULL
#include "NumericsFwd.h"       // for SolverOptions
#include "SOCLCP_Solvers.h"    // for soclcp_VI_ExtraGradient_setDefaultSolv...
#include "SOCLCP_cst.h"        // for SICONOS_SOCLCP_NSGS, SICONOS_SOCLCP_VI_EG
#include "SolverOptions.h"     // for solver_options_nullify, SolverOptions
#include "numerics_verbose.h"  // for numerics_error

int soclcp_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  options->iparam = NULL;
  options->dparam = NULL;
  solver_options_nullify(options);

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
  /* case SICONOS_SOCLCP_NSN_AC: */
  /* { */
  /*   info =    soclcp_AlartCurnier_setDefaultSolverOptions(options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_NSN_FB: */
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
    numerics_error("soclcp_setDefaultSolverOptions", "Unknown Solver");

  }
  }

  return info;
}
