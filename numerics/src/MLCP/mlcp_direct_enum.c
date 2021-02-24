/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*
|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/
#include "mlcp_direct_enum.h"
#ifndef __cplusplus
#include <stdbool.h>                       // for false
#endif
#include <stdio.h>                              // for printf
#include "MLCP_Solvers.h"                       // for mlcp_direct, mlcp_enum
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "SolverOptions.h"                      // for SolverOptions
#include "mlcp_cst.h"                           // for SICONOS_DPARAM_MLCP_S...
#include "mlcp_direct.h"                        // for mlcp_direct_addConfig...
#include "numerics_verbose.h"                   // for numerics_error, verbose


/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"

static int sN;
static int sM;

static int * siWorkEnum = 0;
static int * siWorkDirect = 0;
static double * sdWorkEnum = 0;
static double * sdWorkDirect = 0;

void mlcp_direct_enum_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  sN = problem->n;
  sM = problem->m;
  int iOffset = mlcp_direct_getNbIWork(problem, options);
  int dOffset = mlcp_direct_getNbDWork(problem, options);
  siWorkEnum = options->iWork + iOffset;
  siWorkDirect = options->iWork;
  sdWorkEnum = options->dWork + dOffset;
  sdWorkDirect = options->dWork;
  mlcp_direct_init(problem, options);

}
void mlcp_direct_enum_reset()
{
  mlcp_direct_reset();
  siWorkEnum = 0;
  siWorkDirect = 0;
  sdWorkEnum = 0;
  sdWorkDirect = 0;
}

void mlcp_direct_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("mlcp_direct_enum(...)\n");
  DEBUG_PRINTF("options->iWork = %p\n",  options->iWork);
  if(!siWorkEnum)
  {
    *info = 1;
    numerics_printf_verbose(0,"MLCP_DIRECT_ENUM error, call a non initialised method!!!!!!!!!!!!!!!!!!!!!\n");
    return;
  }
  /*First, try direct solver*/
  options->dWork = sdWorkDirect;
  options->iWork = siWorkDirect;
  mlcp_direct(problem, z, w, info, options);
  if(*info)
  {
    DEBUG_PRINT("Solver direct failed, so run the enum solver\n");
    options->dWork = sdWorkEnum;
    options->iWork = siWorkEnum;
    mlcp_enum(problem, z, w, info, options);
    if(!(*info))
    {
      mlcp_direct_addConfigFromWSolution(problem, w + sN);
    }
    /* Come back to previous memory adress to ensure correct freeing */
    options->dWork = sdWorkDirect;
    options->iWork = siWorkDirect;
  }
  DEBUG_PRINTF("options->iWork = %p\n",  options->iWork);
  DEBUG_END("mlcp_direct_enum(...)\n");
}

void mlcp_direct_enum_set_default(SolverOptions* options)
{
  options->dparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0;
  options->dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12;
  options->dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12;
  options->iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
  options->iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;
  options->filterOn = false;

}
