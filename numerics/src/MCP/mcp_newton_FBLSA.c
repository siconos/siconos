/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include "mcp_newton_FBLSA.h"
#include <assert.h>                       // for assert
#include "FischerBurmeister.h"            // for Jac_F_FB, phi_Mixed_FB
#include "MCP_Solvers.h"                  // for mcp_compute_error, mcp_newt...
#include "MixedComplementarityProblem.h"  // for MixedComplementarityProblem
#include "Newton_methods.h"               // for functions_LSA, init_lsa_fun...
#include "SiconosBlas.h"                  // for cblas_dnrm2
#include "SolverOptions.h"                // for SolverOptions, SICONOS_DPAR...
#include "numerics_verbose.h"             // for numerics_printf

void FB_compute_F_mcp(void* data_opaque, double* z, double* Fmcp)
{
  // Computation of the new value F(z)
  MixedComplementarityProblem* data = (MixedComplementarityProblem *)data_opaque;
  data->compute_Fmcp(data->env, data->n1 + data->n2, z, Fmcp);
}

void FB_compute_H_mcp(void* data_opaque, double* z, double* Fmcp, double* workV1, double* workV2, NumericsMatrix* H)
{
  MixedComplementarityProblem* data = (MixedComplementarityProblem *)data_opaque;

  assert(data->nabla_Fmcp);
  data->compute_nabla_Fmcp(data->env, data->n1 + data->n2, z, data->nabla_Fmcp);

  Jac_F_FB(data->n1, data->n2, z, Fmcp, workV1, workV2, data->nabla_Fmcp, H);
}

void FB_compute_error_mcp(void* data_opaque, double* z, double* w, double* Jac_F_merit, double tol, double* err)
{
  MixedComplementarityProblem* data = (MixedComplementarityProblem *)data_opaque;
  unsigned int n = data->n1 + data->n2;
  err[0] = cblas_dnrm2(n, Jac_F_merit, 1);

  /* If we want to control the convergence of the Newton solver up to the mcp criterion for
     computing error */
  /* mcp_compute_error(data, z , w, err); */

}

void mcp_FB(void* data_opaque, double* z, double* F, double* F_FB)
{
  MixedComplementarityProblem* data = (MixedComplementarityProblem *)data_opaque;
  phi_Mixed_FB(data->n1, data->n2, z, F, F_FB);
}

void mcp_newton_FB_FBLSA(MixedComplementarityProblem* problem, double *z, double* Fmcp, int *info, SolverOptions* options)
{
  numerics_printf("mcp_newton_FB_FBLSA. starts");
  functions_LSA functions_FBLSA_mcp;

  /* This call will set
   * functions_FBLSA_mcp.compute_F to FB_compute_F_mcp
   * functions_FBLSA_mcp.compute_F_merit to mcp_FB
   */
  init_lsa_functions(&functions_FBLSA_mcp, &FB_compute_F_mcp, &mcp_FB);

  /* function to get an element H of T (in our case the clarke Jacobian of F) */
  functions_FBLSA_mcp.compute_H = &FB_compute_H_mcp;

  /* function to compute the error (in our case the norm of the gradient of the merit function) */
  functions_FBLSA_mcp.compute_error = &FB_compute_error_mcp;


  set_lsa_params_data(options, problem->nabla_Fmcp);
  newton_LSA(problem->n1 + problem->n2, z, Fmcp, info, (void *)problem, options, &functions_FBLSA_mcp);

  double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  double  error =0.0;

  mcp_compute_error(problem, z, Fmcp, &error);

  if(error > tolerance)
  {
    numerics_printf("mcp_newton_FB_FBLSA : error = %e > tolerance = %e.", error, tolerance);
    *info = 1;
  }
  else
  {
    numerics_printf("mcp_newton_FB_FBLSA : error = %e < tolerance = %e.", error, tolerance);
    *info = 0;
  }
  options->dparam[SICONOS_DPARAM_RESIDU] = error;

  numerics_printf("mcp_newton_FB_FBLSA. ends");
}

