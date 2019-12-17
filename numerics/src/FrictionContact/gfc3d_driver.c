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
#include <assert.h>                        // for assert
#include <stdio.h>                         // for NULL, fprintf, stderr
#include <stdlib.h>                        // for exit, EXIT_FAILURE
#include "Friction_cst.h"                  // for SICONOS_GLOBAL_FRICTION_3D...
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for gfc3d_driver
#include "NumericsFwd.h"                   // for SolverOptions, GlobalFrict...
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "debug.h"                         // for DEBUG_EXPR
#include "gfc3d_Solvers.h"                 // for gfc3d_ACLMFixedPoint, gfc3...
#include "numerics_verbose.h"              // for numerics_printf_verbose
#include "NumericsMatrix.h"
#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
#endif
int * Global_ipiv = NULL;
int  Global_MisInverse = 0;
int  Global_MisLU = 0;

const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_WR_STR = "GFC3D_NSGS_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR_STR = "GFC3D_NSN_AC_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR_STR = "GFC3D_NSGSV_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_PROX_WR_STR = "GFC3D_PROX_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_DSFP_WR_STR = "GFC3D_DSFP_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_TFP_WR_STR = "GFC3D_TFP_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_STR = "GFC3D_NSGS";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSN_AC_STR = "GFC3D_NSN_AC";
const char* const  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH_STR = "GFC3D_GAMS_PATH";
const char* const  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI_STR = "GFC3D_GAMS_PATHVI";
const char* const  SICONOS_GLOBAL_FRICTION_3D_VI_EG_STR = "GFC3D_VI_EG";
const char* const  SICONOS_GLOBAL_FRICTION_3D_ACLMFP_STR = "GFC3D_ACLMFP";
const char* const  SICONOS_GLOBAL_FRICTION_3D_VI_FPP_STR = "GFC3D_VI_FPP";
const char* const SICONOS_GLOBAL_FRICTION_3D_ADMM_WR_STR = "GFC3D_ADMM_WR";


int gfc3d_driver(GlobalFrictionContactProblem* problem, double *reaction, double *velocity,
                 double* globalVelocity,  SolverOptions* options)
{
  assert(options->isSet);
  DEBUG_EXPR(NV_display(globalVelocity,problem->M->size0););
  if(verbose > 0)
    solver_options_print(options);

  /* Solver name */
  /*  const char* const  name = options->solverName;*/


  int info = -1 ;

  if(problem->dimension != 3)
    numerics_error("gfc3d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* if there is no contact, we compute directly the global velocity as M^{-1}q */
  int m = problem->H->size1;
  if(m ==0)
  {
    numerics_printf_verbose(1,"---- GFC3D - DRIVER . No contact case. Direct computation of global velocity");
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    return 0;
  }

  /* Non Smooth Gauss Seidel (NSGS) */
  switch(options->solverId)
  {
  case SICONOS_GLOBAL_FRICTION_3D_NSGS_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGS_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_nsgs_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGSV_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_nsgs_velocity_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSN_AC_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_nonsmooth_Newton_AlartCurnier_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_PROX_WR:
  {

    numerics_printf_verbose(1," ========================== Call PROX_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_proximal_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_DSFP_WR:
  {

    numerics_printf_verbose(1," ========================== Call DSFP_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_DeSaxceFixedPoint_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_TFP_WR:
  {

    numerics_printf_verbose(1," ========================== Call TFP_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_TrescaFixedPoint_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_NSGS:
  {
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_nsgs(problem, reaction, velocity, globalVelocity,
               &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC:
  {
    gfc3d_nonsmooth_Newton_AlartCurnier(problem, reaction, velocity,
                                        globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH:
  {
    numerics_printf_verbose(1," ========================== Call PATH solver via GAMS for an AVI Friction-Contact 3D problem ==========================\n");
    gfc3d_AVI_gams_path(problem, reaction, velocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI:
  {
    numerics_printf_verbose(1," ========================== Call PATHVI solver via GAMS for an AVI Friction-Contact 3D problem ==========================\n");
    gfc3d_AVI_gams_pathvi(problem, reaction, globalVelocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_VI_FPP:
  {
    gfc3d_VI_FixedPointProjection(problem, reaction, velocity,
                                  globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_VI_EG:
  {
    gfc3d_VI_ExtraGradient(problem, reaction, velocity,
                           globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_ACLMFP:
  {
    gfc3d_ACLMFixedPoint(problem, reaction, velocity,
                         globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_ADMM:
  {
    /* globalFrictionContact_rescaling(problem, 1.0/1.512808e-04, 1.0/1.407230e+01, 1.0); */
    gfc3d_ADMM(problem, reaction, velocity,
               globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_ADMM_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGS_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    Global_ipiv = NULL;
    Global_MisInverse = 0;
    Global_MisLU = 0;
    gfc3d_admm_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  default:
  {
    fprintf(stderr, "Numerics, gfc3d_driver failed. Unknown solver %d.\n", options->solverId);
    exit(EXIT_FAILURE);

  }
  }

  return info;

}

int gfc3d_checkTrivialCaseGlobal(int n, double* q, double* velocity, double* reaction, double * globalVelocity, SolverOptions* options)
{
  /* norm of vector q */
  /*   double qs = cblas_dnrm2( n , q , 1 ); */
  /*   int i; */
  int info = -1;
  /*   if( qs <= DBL_EPSILON )  */
  /*     { */
  /*       // q norm equal to zero (less than DBL_EPSILON) */
  /*       // -> trivial solution: reaction = 0 and velocity = q */
  /*       for( i = 0 ; i < n ; ++i ) */
  /*  { */
  /*    velocity[i] = q[i]; */
  /*    reaction[i] = 0.; */
  /*  } */
  /*       iparam[2] = 0; */
  /*       iparam[4]= 0; */
  /*       dparam[1] = 0.0; */
  /*       dparam[3] = 0.0; */
  /*       info = 0; */
  /*       if(iparam[1]>0) */
  /*  printf("fc3d driver, norm(q) = 0, trivial solution reaction = 0, velocity = q.\n"); */
  /*     } */
  return info;
}
