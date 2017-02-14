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
#include <math.h>
#include <float.h>

#include "MCP_Solvers.h"
#include "SiconosLapack.h"
#include "Newton_methods.h"
#include "FischerBurmeister.h"
#include "min_merit.h"
#include "mcp_newton_FBLSA.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

static void mcp_min(void* data_opaque, double* z, double* F, double* Fmin)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;

  F_min(data->n1, data->n2, z, F, Fmin);
}

static void min_compute_H_mcp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, NumericsMatrix* H)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;

  data->compute_nabla_Fmcp(data->env, data->n1 + data->n2, z, data->nabla_Fmcp);

  Jac_F_min(data->n1, data->n2, z, F, data->nabla_Fmcp, H);
}

void mcp_newton_minFBLSA(MixedComplementarityProblem2* problem, double *z, double* Fmcp, int *info , SolverOptions* options)
{
  functions_LSA functions_minFBLSA_mcp;
  init_lsa_functions(&functions_minFBLSA_mcp, &FB_compute_F_mcp, &mcp_FB);
  functions_minFBLSA_mcp.compute_H = &FB_compute_H_mcp;
  functions_minFBLSA_mcp.compute_error = &FB_compute_error_mcp;
  functions_minFBLSA_mcp.compute_RHS_desc = &mcp_min;
  functions_minFBLSA_mcp.compute_H_desc = &min_compute_H_mcp;

  set_lsa_params_data(options, problem->nabla_Fmcp);
  newton_LSA(problem->n1 + problem->n2, z, Fmcp, info, (void *)problem, options, &functions_minFBLSA_mcp);
}
