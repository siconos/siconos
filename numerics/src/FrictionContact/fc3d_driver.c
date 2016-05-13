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

#include "NumericsOptions.h"
#include "fc3d_Solvers.h"
#include "NonSmoothDrivers.h"

char *  SICONOS_FRICTION_3D_NSGS_STR = "FC3D_NSGS";
char *  SICONOS_FRICTION_3D_NSGSV_STR = "FC3D_NSGSV";
char *  SICONOS_FRICTION_3D_TFP_STR = "FC3D_TFP";
char *  SICONOS_FRICTION_3D_LOCALAC_STR = "FC3D_LOCALAC";
char *  SICONOS_FRICTION_3D_LOCALFB_STR = "FC3D_LOCALFB";
char *  SICONOS_FRICTION_3D_NSN_AC_STR = "FC3D_NSN_AC";
char *  SICONOS_FRICTION_3D_NSN_FB_STR = "FC3D_NSN_FB";
char *  SICONOS_FRICTION_3D_NSN_NM_STR = "FC3D_NSN_NM";
char *  SICONOS_FRICTION_3D_DSFP_STR = "FC3D_DeSaxceFixedPoint";
char *  SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint_STR = "FC3D_NCPGlockerFBFixedPoint";
char *  SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_STR = "FC3D_ONECONTACT_NSN_AC";
char *  SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP_STR = "FC3D_ONECONTACT_NSN_AC_GP";
char *  SICONOS_FRICTION_3D_NCPGlockerFBNewton_STR = "FC3D_NCPGlockerFBNewton";
char * SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization_STR = "FC3D_ProjectionOnConeWithDiagonalization";
char * SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_STR = "FC3D_ProjectionOnCone";
char * SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration_STR = "FC3D_ProjectionOnConeWithLocalIteration";
char * SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization_STR = "FC3D_projectionOnConeWithRegularization";
char * SICONOS_FRICTION_3D_NCPGlockerFBPATH_STR = "FC3D_NCPGlockerFBPATH";
char * SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder_STR = "FC3D_projectionOnCylinder";
char * SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration_STR =  "FC3D_projectionOnCylinderWithLocalIteration";
char * SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity_STR = "FC3D_ProjectionOnCone_velocity";
char * SICONOS_FRICTION_3D_PGoC_STR = "FC3D_PGoC";
char * SICONOS_FRICTION_3D_DeSaxceFixedPoint_STR = "FC3D_DeSaxceFixedPoint";
char * SICONOS_FRICTION_3D_EG_STR = "FC3D_ExtraGradient";
char * SICONOS_FRICTION_3D_FPP_STR = "FC3D_FixedPointProjection";
char * SICONOS_FRICTION_3D_VI_EG_STR = "FC3D_VI_ExtraGradient";
char * SICONOS_FRICTION_3D_VI_FPP_STR = "FC3D_VI_FixedPointProjection";
char * SICONOS_FRICTION_3D_HP_STR = "FC3D_HyperplaneProjection";
char * SICONOS_FRICTION_3D_PROX_STR = "FC3D_PROX";
char * SICONOS_FRICTION_3D_GAMS_PATH_STR = "FC3D_GAMS_PATH";
char * SICONOS_FRICTION_3D_GAMS_PATHVI_STR = "FC3D_GAMS_PATHVI";
char * SICONOS_FRICTION_3D_GAMS_LCP_PATH_STR = "FC3D_GAMS_LCP_PATH";
char * SICONOS_FRICTION_3D_GAMS_LCP_PATHVI_STR = "FC3D_GAMS_LCP_PATHVI";
char * SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_STR = "FC3D_QUARTIC";
char * SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU_STR = "FC3D_QUARTIC_NU";
char *  SICONOS_FRICTION_3D_ACLMFP_STR = "FC3D_ACLMFP";
char *  SICONOS_FRICTION_3D_SOCLCP_STR = "FC3D_SOCLCP";


void snPrintf(int level, SolverOptions* opts, const char *fmt, ...);

int fc3d_driver(FrictionContactProblem* problem,
                             double *reaction, double *velocity,
                             SolverOptions* options,
                             NumericsOptions* global_options)
{
  if (options == NULL)
    numericsError("fc3d_driver", "null input for solver and/or global options");

  int setnumericsoptions=0;

  /* Set global options */
  if (global_options)
  {
    setNumericsOptions(global_options);
    options->numericsOptions = (NumericsOptions*) malloc(sizeof(NumericsOptions));
    options->numericsOptions->verboseMode = global_options->verboseMode;
    setnumericsoptions=1;
  }

  int NoDefaultOptions = options->isSet; /* true(1) if the SolverOptions structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(3, options);

  if (verbose > 1)
    printSolverOptions(options);

  /* Solver name */
  /*char * name = options->solverName;*/

  int info = -1 ;

  if (problem->dimension != 3)
    numericsError("fc3d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* Check for trivial case */
  info = checkTrivialCase(problem, velocity, reaction, options);


  if (info == 0)
    goto exit;


  switch (options->solverId)
  {
    /* Non Smooth Gauss Seidel (NSGS) */
  case SICONOS_FRICTION_3D_NSGS:
  {
    snPrintf(1, options,
             " ========================== Call NSGS solver for Friction-Contact 3D problem ==========================\n");
    fc3d_nsgs(problem, reaction , velocity , &info , options);
    break;
  }
  case SICONOS_FRICTION_3D_NSGSV:
  {
    snPrintf(1, options,
             " ========================== Call NSGSV solver for Friction-Contact 3D problem ==========================\n");
    fc3d_nsgs_velocity(problem, reaction , velocity , &info , options);
    break;
  }
  /* Proximal point algorithm */
  case SICONOS_FRICTION_3D_PROX:
  {
    snPrintf(1, options,
             " ========================== Call PROX (Proximal Point) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_proximal(problem, reaction , velocity , &info , options);
    break;
  }
  /* Tresca Fixed point algorithm */
  case SICONOS_FRICTION_3D_TFP:
  {
    snPrintf(1, options,
             " ========================== Call TFP (Tresca Fixed Point) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_TrescaFixedPoint(problem, reaction , velocity , &info , options);
    break;
  }
  /* ACLM Fixed point algorithm */
  case SICONOS_FRICTION_3D_ACLMFP:
  {
    snPrintf(1, options,
             " ========================== Call ACLM (Acary Cadoux Lemarechal Malick Fixed Point) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_ACLMFixedPoint(problem, reaction , velocity , &info , options);
    break;
  }
  /* SOCLCP Fixed point algorithm */
  case SICONOS_FRICTION_3D_SOCLCP:
  {
    snPrintf(1, options,
             " ========================== Call SOCLCP solver for Friction-Contact 3D problem (Associated one) ==========================\n");
    fc3d_SOCLCP(problem, reaction , velocity , &info , options);
    break;
  }
  /* De Saxce Fixed point algorithm */
  case SICONOS_FRICTION_3D_DSFP:
  {
    snPrintf(1, options,
            " ========================== Call DeSaxce Fixed Point (DSFP) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_DeSaxceFixedPoint(problem, reaction , velocity , &info , options);
    break;
  }
  /* Fixed point projection algorithm */
  case SICONOS_FRICTION_3D_FPP:
  {
    snPrintf(1, options,
            " ========================== Call Fixed Point Projection (FPP) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_fixedPointProjection(problem, reaction , velocity , &info , options);
    break;
  }

  /* Extra Gradient algorithm */
  case SICONOS_FRICTION_3D_EG:
  {
    snPrintf(1, options,
            " ========================== Call ExtraGradient (EG) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_ExtraGradient(problem, reaction , velocity , &info , options);
    break;
  }
  /* VI Fixed Point Projection algorithm */
  case SICONOS_FRICTION_3D_VI_FPP:
  {
    snPrintf(1, options,
            " ========================== Call VI_FixedPointProjection (VI_FPP) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_VI_FixedPointProjection(problem, reaction , velocity , &info , options);
    break;
  }
  /* VI Extra Gradient algorithm */
  case SICONOS_FRICTION_3D_VI_EG:
  {
    snPrintf(1, options,
            " ========================== Call VI_ExtraGradient (VI_EG) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_VI_ExtraGradient(problem, reaction , velocity , &info , options);
    break;
  }
  /* Hyperplane Projection algorithm */
  case SICONOS_FRICTION_3D_HP:
  {
    snPrintf(1, options,
            " ========================== Call Hyperplane Projection (HP) solver for Friction-Contact 3D problem ==========================\n");
    fc3d_HyperplaneProjection(problem, reaction , velocity , &info , options);
    break;
  }
  /* Alart Curnier in local coordinates */
  case SICONOS_FRICTION_3D_NSN_AC:
  {
    snPrintf(1, options,
            " ========================== Call Alart Curnier solver for Friction-Contact 3D problem ==========================\n");
    if (problem->M->matrix0)
    {
      fc3d_nonsmooth_Newton_AlartCurnier(problem, reaction , velocity , &info , options);
    }
    else
    {
      fc3d_nonsmooth_Newton_AlartCurnier(problem, reaction , velocity , &info , options);
    }
    break;
  }
  /* Fischer Burmeister in local coordinates */
  case SICONOS_FRICTION_3D_NSN_FB:
  {
    snPrintf(1, options,
            " ========================== Call Fischer Burmeister solver for Friction-Contact 3D problem ==========================\n");
    fc3d_nonsmooth_Newton_FischerBurmeister(problem, reaction , velocity , &info , options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_NM:
  {
    snPrintf(1, options,
            " ========================== Call natural map solver for Friction-Contact 3D problem ==========================\n");
    fc3d_nonsmooth_Newton_NaturalMap(problem, reaction , velocity , &info , options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
  {
    snPrintf(1, options,
            " ========================== Call Quartic solver for Friction-Contact 3D problem ==========================\n");
    fc3d_unitary_enumerative(problem, reaction , velocity , &info , options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC:
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP:
  {
    snPrintf(1, options,
            " ========================== Call Newton-based solver for one contact Friction-Contact 3D problem ==========================\n");
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, problem, options);
    info = fc3d_onecontact_nonsmooth_Newton_solvers_solve(problem, reaction , options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  {
    snPrintf(1, options,
            " ========================== Call Projection on cone solver for one contact Friction-Contact 3D problem ==========================\n");
    fc3d_projectionOnConeWithLocalIteration_initialize(problem, problem, options);
    info = fc3d_projectionOnConeWithLocalIteration_solve(problem, reaction , options);
    fc3d_projectionOnConeWithLocalIteration_free(problem, problem, options);

    break;
  }
  case SICONOS_FRICTION_3D_GAMS_PATH:
  {
    snPrintf(1, options,
            " ========================== Call PATH solver via GAMS for an AVI Friction-Contact 3D problem ==========================\n");
    fc3d_AVI_gams_path(problem, reaction , velocity, &info, options);
    break;
  }
  case SICONOS_FRICTION_3D_GAMS_PATHVI:
  {
    snPrintf(1, options,
            " ========================== Call PATHVI solver via GAMS for an AVI Friction-Contact 3D problem ==========================\n");
    fc3d_AVI_gams_pathvi(problem, reaction , velocity, &info, options);
    break;
  }
  case SICONOS_FRICTION_3D_GAMS_LCP_PATH:
  {
    snPrintf(1, options,
            " ========================== Call PATH solver via GAMS for an LCP-based reformulation of the AVI Friction-Contact 3D problem ==========================\n");
    fc3d_lcp_gams_path(problem, reaction , velocity, &info, options);
    break;
  }
  case SICONOS_FRICTION_3D_GAMS_LCP_PATHVI:
  {
    snPrintf(1, options,
            " ========================== Call PATHVI solver via GAMS for an LCP-based reformulation of the AVI Friction-Contact 3D problem ==========================\n");
    fc3d_lcp_gams_pathvi(problem, reaction , velocity, &info, options);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, fc3d_driver failed. Unknown solver.\n");
    exit(EXIT_FAILURE);

  }
  }

exit:

  if (setnumericsoptions)
  {
      free(options->numericsOptions);
      options->numericsOptions = NULL;
  }

  return info;

}

int checkTrivialCase(FrictionContactProblem* problem, double* velocity,
                     double* reaction, SolverOptions* options)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  /* Dimension of the problem */
  int n = 3 * nc;
  int i = 0;
  /*take off? R=0 ?*/
  for (i = 0; i < nc; i++)
  {
    if (q[3 * i] < -DBL_EPSILON)
      return -1;
  }
  for (i = 0 ; i < n ; ++i)
  {
    velocity[i] = q[i];
    reaction[i] = 0.;
  }
  options->iparam[2] = 0;
  options->iparam[4] = 0;
  options->dparam[1] = 0.0;
  options->dparam[3] = 0.0;
  if (verbose == 1)
    printf("fc3d driver, take off, trivial solution reaction = 0, velocity = q.\n");
  return 0;
}

#include <stdarg.h>
void snPrintf(int level, SolverOptions* opts, const char *fmt, ...)
{
  assert(opts);
  assert(opts->numericsOptions);

  NumericsOptions* nopts = opts->numericsOptions;

  if (nopts && nopts->verboseMode >= level)
  {
    va_list args;
    va_start(args,fmt);
    printf("Siconos/Numerics: ");
    vprintf(fmt,args);
    va_end(args);
  }
}
