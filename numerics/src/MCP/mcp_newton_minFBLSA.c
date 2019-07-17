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


#include <stdio.h>
#include <math.h>
#include <float.h>

#include "MCP_Solvers.h"
#include "MCP_cst.h"
#include "SiconosLapack.h"
#include "Newton_methods.h"
#include "FischerBurmeister.h"
#include "min_merit.h"
#include "numerics_verbose.h"
#include "mcp_newton_FBLSA.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

static void mcp_min(void* data_opaque, double* z, double* F, double* Fmin)
{
  MixedComplementarityProblem* data = (MixedComplementarityProblem *)data_opaque;

  F_min(data->n1, data->n2, z, F, Fmin);
}

static void min_compute_H_mcp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, NumericsMatrix* H)
{
  MixedComplementarityProblem* data = (MixedComplementarityProblem *)data_opaque;

  data->compute_nabla_Fmcp(data->env, data->n1 + data->n2, z, data->nabla_Fmcp);

  Jac_F_min(data->n1, data->n2, z, F, data->nabla_Fmcp, H);
}

void mcp_newton_min_FBLSA(MixedComplementarityProblem* problem, double *z, double* Fmcp, int *info , SolverOptions* options)
{
  numerics_printf("mcp_newton_min_FBLSA. starts");
  functions_LSA functions_minFBLSA_mcp;
  init_lsa_functions(&functions_minFBLSA_mcp, &FB_compute_F_mcp, &mcp_FB);
  functions_minFBLSA_mcp.compute_H = &FB_compute_H_mcp;
  functions_minFBLSA_mcp.compute_error = &FB_compute_error_mcp;
  functions_minFBLSA_mcp.compute_RHS_desc = &mcp_min;
  functions_minFBLSA_mcp.compute_H_desc = &min_compute_H_mcp;
  
  options->internalSolvers->dparam[0] = options->dparam[0];
  options->internalSolvers->iparam[0] = options->iparam[0];

  set_lsa_params_data(options->internalSolvers, problem->nabla_Fmcp);
  newton_LSA(problem->n1 + problem->n2, z, Fmcp, info, (void *)problem,
             options->internalSolvers,
             &functions_minFBLSA_mcp);
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  double  error =0.0;

  mcp_compute_error(problem, z , Fmcp, &error);

  if (error > tolerance)
  {
    numerics_printf("mcp_newton_min_FBLSA : error = %e > tolerance = %e.", error, tolerance);
    *info = 1;
  }
  else
  {
    numerics_printf("mcp_newton_min_FBLSA : error = %e < tolerance = %e.", error, tolerance);
    *info = 0;
  }
 
  
  options->iparam[SICONOS_IPARAM_ITER_DONE] = options->internalSolvers->iparam[SICONOS_IPARAM_ITER_DONE];
  options->dparam[SICONOS_DPARAM_RESIDU] = error;

  numerics_printf("mcp_newton_min_FBLSA. ends");
}

int mcp_newton_min_FBLSA_setDefaultSolverOptions(
  MixedComplementarityProblem* problem,
  SolverOptions* options)
{
  numerics_printf_verbose(1,"mcp_newton_min_FBLSA_setDefaultSolverOptions");

  options->solverId = SICONOS_MCP_NEWTON_MIN_FBLSA;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  
  newton_lsa_setDefaultSolverOptions(options->internalSolvers);

  options->internalSolvers->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->internalSolvers->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  return 0;
}
