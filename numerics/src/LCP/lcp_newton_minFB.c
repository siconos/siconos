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

#include <assert.h>  // for assert

#include "LCP_Solvers.h"                   // for lcp_newton_minFB
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NonSmoothSolvers/Newton_methods.h"                // for functions_LSA, init_lsa_fu...
#include "NumericsFwd.h"                   // for LinearComplementarityProblem
#include "lcp_newton_FB.h"                 // for FB_compute_F_lcp, FB_compu...
#include "lcp_cst.h"                       // for SICONOS_LCP_NEWTON_MIN_FBLSA
#include "min_merit.h"                     // for F_min, Jac_F_min

/* Solver registration system */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"
static void lcp_min(void* data_opaque, double* z, double* F, double* Fmin) {
  F_min(0, ((LinearComplementarityProblem*)data_opaque)->size, z, F, Fmin);
}

static void min_compute_H_lcp(void* data_opaque, double* z, double* F, double* workV1,
                              double* workV2, NumericsMatrix* H) {
  LinearComplementarityProblem* data = (LinearComplementarityProblem*)data_opaque;
  unsigned int n = data->size;
  assert(data->M);

  Jac_F_min(0, n, z, F, data->M, H);
}

void lcp_newton_minFB(LinearComplementarityProblem* problem, double* z, double* w, int* info,
                      SolverOptions* options) {
  functions_LSA functions_minFBLSA_lcp;
  init_lsa_functions(&functions_minFBLSA_lcp, &FB_compute_F_lcp, &lcp_FB);
  functions_minFBLSA_lcp.compute_H = &FB_compute_H_lcp;
  functions_minFBLSA_lcp.compute_error = &FB_compute_error_lcp;
  functions_minFBLSA_lcp.compute_RHS_desc = &lcp_min;
  functions_minFBLSA_lcp.compute_H_desc = &min_compute_H_lcp;

  set_lsa_params_data(options, problem->M);
  newton_LSA(problem->size, z, w, info, (void*)problem, options, &functions_minFBLSA_lcp);
}

static void lcp_newton_minFB_set_default(SolverOptions* options) {
  /* No specific defaults needed */
  (void)options;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_LCP_NEWTON_MIN_FBLSA in the global solver registry.
 */

static int lcp_newton_minFB_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  lcp_newton_minFB_set_default(options);
  return NUMERICS_OK;
}

static int lcp_newton_minFB_solve_wrap(void* problem, double* reaction,
                                       double* velocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  lcp_newton_minFB((LinearComplementarityProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void lcp_newton_minFB_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_LCP_NEWTON_MIN_FBLSA, "LCP_NEWTON_MIN_FBLSA",
                       "Newton Min FBLSA solver for LCP",
                       lcp_newton_minFB_init_wrap,
                       lcp_newton_minFB_solve_wrap,
                       lcp_newton_minFB_free_wrap,
                       NULL,  /* error function */
                       lcp_newton_minFB_set_default,  /* set_default */
                       1000,  /* default_max_iter */
                       1e-6,  /* default_tol */
                       0);     /* is_local_solver */
