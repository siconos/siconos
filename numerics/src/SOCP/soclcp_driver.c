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
#include <assert.h>                                       // for assert
#include <float.h>                                        // for DBL_EPSILON
#include <stdio.h>                                        // for fprintf
#include <stdlib.h>                                       // for exit, EXIT_...
#include "NonSmoothDrivers.h"                             // for soclcp_driver
#include "NumericsFwd.h"                                  // for SolverOptions
#include "SOCLCP_Solvers.h"                               // for soclcp_VI_E...
#include "SOCLCP_cst.h"                                   // for SICONOS_SOC...
#include "SecondOrderConeLinearComplementarityProblem.h"  // for SecondOrder...
#include "SolverOptions.h"                                // for SolverOptions
#include "numerics_verbose.h"                             // for numerics_pr...

const char* const   SICONOS_SOCLCP_NSGS_STR = "SOCLCP_NSGS";
const char* const   SICONOS_SOCLCP_NSGSV_STR = "SOCLCP_NSGSV";
const char* const   SICONOS_SOCLCP_TFP_STR = "SOCLCP_TFP";
const char* const   SICONOS_SOCLCP_NSN_AC_STR = "SOCLCP_NSN_AC";
const char* const   SICONOS_SOCLCP_NSN_FB_STR = "SOCLCP_NSN_FB";
const char* const   SICONOS_SOCLCP_DSFP_STR = "SOCLCP_DeSaxceFixedPoint";
const char* const   SICONOS_SOCLCP_NCPGlockerFBFixedPoint_STR = "SOCLCP_NCPGlockerFBFixedPoint";
const char* const   SICONOS_SOCLCP_AlartCurnierNewton_STR = "SOCLCP_AlartCurnierNewton";
const char* const   SICONOS_SOCLCP_DampedAlartCurnierNewton_STR = "SOCLCP_DampedAlartCurnierNewton";
const char* const   SICONOS_SOCLCP_NCPGlockerFBNewton_STR = "SOCLCP_NCPGlockerFBNewton";
const char* const  SICONOS_SOCLCP_ProjectionOnConeWithDiagonalization_STR = "SOCLCP_ProjectionOnConeWithDiagonalization";
const char* const  SICONOS_SOCLCP_ProjectionOnCone_STR = "SOCLCP_ProjectionOnCone";
const char* const  SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration_STR = "SOCLCP_ProjectionOnConeWithLocalIteration";
const char* const  SICONOS_SOCLCP_ProjectionOnConeWithRegularization_STR = "SOCLCP_ProjectionOnConeWithRegularization";
const char* const  SICONOS_SOCLCP_NCPGlockerFBPATH_STR = "SOCLCP_NCPGlockerFBPATH";
const char* const  SICONOS_SOCLCP_projectionOnCylinder_STR = "SOCLCP_projectionOnCylinder";
const char* const  SICONOS_SOCLCP_ProjectionOnCone_velocity_STR = "SOCLCP_ProjectionOnCone_velocity";
const char* const  SICONOS_SOCLCP_PGoC_STR = "SOCLCP_PGoC";
const char* const  SICONOS_SOCLCP_DeSaxceFixedPoint_STR = "SOCLCP_DeSaxceFixedPoint";
const char* const  SICONOS_SOCLCP_EG_STR = "SOCLCP_ExtraGradient";
const char* const  SICONOS_SOCLCP_FPP_STR = "SOCLCP_FixedPointProjection";
const char* const  SICONOS_SOCLCP_VI_EG_STR = "SOCLCP_VI_ExtraGradient";
const char* const  SICONOS_SOCLCP_VI_FPP_STR = "SOCLCP_VI_FixedPointProjection";
const char* const  SICONOS_SOCLCP_HP_STR = "SOCLCP_HyperplaneProjection";
const char* const  SICONOS_SOCLCP_PROX_STR = "SOCLCP_PROX";
const char* const  SICONOS_SOCLCP_QUARTIC_STR = "SOCLCP_QUARTIC";
const char* const  SICONOS_SOCLCP_QUARTIC_NU_STR = "SOCLCP_QUARTIC_NU";


int soclcp_driver(SecondOrderConeLinearComplementarityProblem* problem,
                  double *r, double *v,
                  SolverOptions* options)
{
  if(options == NULL)
    numerics_error("soclcp_driver", "null input for solver and/or global options");

  assert(options->isSet);

  if(verbose > 0)
    solver_options_print(options);

  /* Solver name */
  /*const char* const  name = options->solverName;*/

  int info = -1 ;

  /* Check for trivial case */
  info = soclcp_checkTrivialCase(problem, v, r, options);


  if(info == 0)
    return info;


  switch(options->solverId)
  {
  /* Non Smooth Gauss Seidel (NSGS) */
  case SICONOS_SOCLCP_NSGS:
  {
    numerics_printf(" ========================== Call NSGS solver for Second Order Cone LCP problem ==========================\n");
    soclcp_nsgs(problem, r , v , &info , options);
    break;
  }
  /* case SICONOS_SOCLCP_NSGSV: */
  /* { */
  /*   numerics_printf(numerics_printf(" ========================== Call NSGSV solver for Second Order Cone LCP problem ==========================\n"); */
  /*   soclcp_nsgs_v(problem, r , v , &info , options); */
  /*   break; */
  /* } */
  /* /\* Proximal point algorithm *\/ */
  /* case SICONOS_SOCLCP_PROX: */
  /* { */
  /*   numerics_printf(numerics_printf(" ========================== Call PROX (Proximal Point) solver for Second Order Cone LCP problem ==========================\n"); */
  /*   soclcp_proximal(problem, r , v , &info , options); */
  /*   break; */
  /* } */
  /* /\* Tresca Fixed point algorithm *\/ */
  /* case SICONOS_SOCLCP_TFP: */
  /* { */
  /*   numerics_printf(" ========================== Call TFP (Tresca Fixed Point) solver for Second Order Cone LCP problem ==========================\n"); */
  /*   soclcp_TrescaFixedPoint(problem, r , v , &info , options); */
  /*   break; */
  /* } */
  case SICONOS_SOCLCP_VI_FPP:
  {
    numerics_printf(" ========================== Call VI_FixedPointProjection (VI_FPP) solver for Second Order Cone LCP problem ==========================\n");
    soclcp_VI_FixedPointProjection(problem, r , v , &info , options);
    break;
  }
  /* VI Extra Gradient algorithm */
  case SICONOS_SOCLCP_VI_EG:
  {
    numerics_printf(" ========================== Call VI_ExtraGradient (VI_EG) solver for Second Order Cone LCP problem ==========================\n");
    soclcp_VI_ExtraGradient(problem, r , v , &info , options);
    break;
  }
  /* /\* Hyperplane Projection algorithm *\/ */
  /* case SICONOS_SOCLCP_HP: */
  /* { */
  /*   numerics_printf(" ========================== Call Hyperplane Projection (HP) solver for Second Order Cone LCP problem ==========================\n"); */
  /*   soclcp_HyperplaneProjection(problem, r , v , &info , options); */
  /*   break; */
  /* } */
  /* /\* Alart Curnier in local coordinates *\/ */
  /* case SICONOS_SOCLCP_NSN_AC: */
  /* { */
  /*   numerics_printf(" ========================== Call Alart Curnier solver for Second Order Cone LCP problem ==========================\n"); */
  /*   if (problem->M->matrix0) */
  /*   { */
  /*     soclcp_nonsmooth_Newton_AlartCurnier(problem, r , v , &info , options); */
  /*   } */
  /*   else */
  /*   { */
  /*     soclcp_nonsmooth_Newton_AlartCurnier(problem, r , v , &info , options); */
  /*   } */
  /*   break; */
  /* } */
  /* /\* Fischer Burmeister in local coordinates *\/ */
  /* case SICONOS_SOCLCP_NSN_FB: */
  /* { */
  /*   numerics_printf(" ========================== Call Fischer Burmeister solver for Second Order Cone LCP problem ==========================\n"); */
  /*   soclcp_nonsmooth_Newton_FischerBurmeister(problem, r , v , &info , options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_QUARTIC_NU: */
  /* case SICONOS_SOCLCP_QUARTIC: */
  /* { */
  /*   numerics_printf(" ========================== Call Quartic solver for Second Order Cone LCP problem ==========================\n"); */
  /*   soclcp_unitary_enumerative(problem, r , v , &info , options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_AlartCurnierNewton: */
  /* case SICONOS_SOCLCP_DampedAlartCurnierNewton: */
  /* { */
  /*   numerics_printf(" ========================== Call Quartic solver for Second Order Cone LCP problem ==========================\n"); */
  /*   info =soclcp_Newton_solve(problem, r , options); */
  /*   break; */
  /* } */
  default:
  {
    fprintf(stderr, "Numerics, SecondOrderConeLinearComplementarity_driver failed. Unknown solver.\n");
    exit(EXIT_FAILURE);

  }
  }

  return info;

}

int soclcp_checkTrivialCase(SecondOrderConeLinearComplementarityProblem* problem, double* v,
                            double* r, SolverOptions* options)
{
  /* Number of contacts */
  int nc = problem->nc;
  double* q = problem->q;

  int i = 0;
  /* Dimension of the problem */
  int n = problem->n;

  /* r=0 ?*/
  int normal=0;
  for(i = 0; i < nc; i++)
  {
    normal = problem->coneIndex[i];
    if(q[normal] < -DBL_EPSILON)
      return -1;

  }
  for(i = 0 ; i < n ; ++i)
  {
    v[i] = q[i];
    r[i] = 0.;
  }
  options->iparam[2] = 0;
  options->iparam[4] = 0;
  options->dparam[SICONOS_IPARAM_ITER_DONE] = 0.0;
  options->dparam[3] = 0.0;
  if(verbose == 1)
    printf("SecondOrderConeLinearComplementarity driver, trivial solution r = 0, v = q.\n");
  return 0;
}
