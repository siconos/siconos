/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "MCP_Solvers.h"                  // for mcp_compute_error, mcp_newt...
#include "MCP_cst.h"
#include "MixedComplementarityProblem.h"  // for MixedComplementarityProblem
#include "Newton_methods.h"               // for functions_LSA, init_lsa_fun...
#include "NumericsFwd.h"                  // for MixedComplementarityProblem
#include "SolverOptions.h"                // for SolverOptions, SICONOS_DPAR...
#include "mcp_newton_FBLSA.h"             // for FB_compute_F_mcp, FB_comput...
#include "min_merit.h"                    // for F_min, Jac_F_min
#include "numerics_verbose.h"

/* Solver registration system */
#include "solver_registry.h"
#include "numerics_errors.h"

static void mcp_min(void* data_opaque, double* z, double* F, double* Fmin) {
  MixedComplementarityProblem* data = (MixedComplementarityProblem*)data_opaque;

  F_min(data->n1, data->n2, z, F, Fmin);
}

static void min_compute_H_mcp(void* data_opaque, double* z, double* F, double* workV1,
                              double* workV2, NumericsMatrix* H) {
  MixedComplementarityProblem* data = (MixedComplementarityProblem*)data_opaque;

  data->compute_nabla_Fmcp(data->env, data->n1 + data->n2, z, data->nabla_Fmcp);

  Jac_F_min(data->n1, data->n2, z, F, data->nabla_Fmcp, H);
}

void mcp_newton_min_FBLSA(MixedComplementarityProblem* problem, double* z, double* Fmcp,
                          int* info, SolverOptions* options) {
  numerics_printf("mcp_newton_min_FBLSA. starts");
  functions_LSA functions_minFBLSA_mcp;
  init_lsa_functions(&functions_minFBLSA_mcp, &FB_compute_F_mcp, &mcp_FB);
  functions_minFBLSA_mcp.compute_H = &FB_compute_H_mcp;
  functions_minFBLSA_mcp.compute_error = &FB_compute_error_mcp;
  functions_minFBLSA_mcp.compute_RHS_desc = &mcp_min;
  functions_minFBLSA_mcp.compute_H_desc = &min_compute_H_mcp;

  set_lsa_params_data(options, problem->nabla_Fmcp);
  newton_LSA(problem->n1 + problem->n2, z, Fmcp, info, (void*)problem, options,
             &functions_minFBLSA_mcp);
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  double error = 0.0;

  mcp_compute_error(problem, z, Fmcp, &error);

  if (error > tolerance) {
    numerics_printf("mcp_newton_min_FBLSA : error = %e > tolerance = %e.", error, tolerance);
    *info = 1;
  } else {
    numerics_printf("mcp_newton_min_FBLSA : error = %e < tolerance = %e.", error, tolerance);
    *info = 0;
  }

  options->dparam[SICONOS_DPARAM_RESIDU] = error;

  numerics_printf("mcp_newton_min_FBLSA. ends");
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers MCP_NEWTON_MIN_FBLSA in the global solver registry.
 */

static void mcp_newton_minfblsa_set_default(SolverOptions* options) {
  SOLVER_MAX_ITER(options) = 1000;
  SOLVER_TOL(options) = 1e-4;
}

static int mcp_newton_minfblsa_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int mcp_newton_minfblsa_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  MixedComplementarityProblem* mcp = (MixedComplementarityProblem*)problem;
  mcp_newton_min_FBLSA(mcp, z, w, &info, options);
  return info;
}

static void mcp_newton_minfblsa_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_MCP_NEWTON_MIN_FBLSA, "MCP_NEWTON_MIN_FBLSA",
                       "Newton minFBLSA solver for Mixed Complementarity Problems",
                       mcp_newton_minfblsa_init_wrap, mcp_newton_minfblsa_solve_wrap, mcp_newton_minfblsa_free_wrap, NULL,
                       mcp_newton_minfblsa_set_default, 1000, 1e-4, 0);
