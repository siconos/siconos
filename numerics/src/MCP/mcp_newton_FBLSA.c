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

#include "MCP_Solvers.h"
#include "MCP_cst.h"
#include "SiconosBlas.h"
#include "Newton_Methods.h"
#include "FischerBurmeister.h"

#include "mcp_newton_FBLSA.h"

void FB_compute_F_mcp(void* data_opaque, double* z, double* Fmcp)
{
  // Computation of the new value F(z)
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;
  data->compute_Fmcp(data->env, data->n1, data->n2, z, Fmcp);
}

void FB_compute_H_mcp(void* data_opaque, double* z, double* Fmcp, double* workV1, double* workV2, NumericsMatrix* H)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;

  assert(data->nabla_Fmcp);
  data->compute_nabla_Fmcp(data->env, data->n1, data->n2, z, data->nabla_Fmcp);

  Jac_F_FB(data->n1, data->n2, z, Fmcp, workV1, workV2, data->nabla_Fmcp, H);
}

void FB_compute_error_mcp(void* data_opaque, double* z, double* w, double* Jac_F_merit, double tol, double* err)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;
  unsigned int n = data->n1 + data->n2;
  err[0] = cblas_dnrm2(n, Jac_F_merit, 1);
}

void mcp_FB(void* data_opaque, double* z, double* F, double* F_FB)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;
  phi_Mixed_FB(data->n1, data->n2, z, F, F_FB);
}

void mcp_newton_FBLSA(MixedComplementarityProblem2* problem, double *z, double* Fmcp, int *info , SolverOptions* options)
{
  functions_LSA functions_FBLSA_mcp;
  init_lsa_functions(&functions_FBLSA_mcp, &FB_compute_F_mcp, &mcp_FB);
  functions_FBLSA_mcp.compute_H = &FB_compute_H_mcp;
  functions_FBLSA_mcp.compute_error = &FB_compute_error_mcp;

  set_lsa_params_data(options, problem->nabla_Fmcp);
  newton_LSA(problem->n1 + problem->n2, z, Fmcp, info, (void *)problem, options, &functions_FBLSA_mcp);
}
