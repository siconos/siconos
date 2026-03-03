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

#include "ncp_newton_FBLSA.h"  // for FB_compute_F_ncp, FB_co...

#include "FischerBurmeister.h"                // for Jac_F_FB, phi_FB
#include "NCP_Solvers.h"                      // for ncp_newton_FBLSA
#include "NCP_cst.h"
#include "NonSmoothSolvers/Newton_methods.h"                   // for functions_LSA, init_lsa...
#include "NonlinearComplementarityProblem.h"  // for NonlinearComplementarit...
#include "NumericsFwd.h"                      // for NonlinearComplementarit...
#include "SiconosBlas.h"                      // for cblas_dnrm2

/* Solver registration system */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"
void ncp_FB(void* data_opaque, double* z, double* F, double* F_FB) {
  phi_FB(((NonlinearComplementarityProblem*)data_opaque)->n, z, F, F_FB);
}

void FB_compute_F_ncp(void* data_opaque, double* z, double* F) {
  // Computation of the new value F(z)
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem*)data_opaque;
  data->compute_F(data->env, data->n, z, F);
}

void FB_compute_H_ncp(void* data_opaque, double* z, double* F, double* workV1, double* workV2,
                      NumericsMatrix* H) {
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem*)data_opaque;

  data->compute_nabla_F(data->env, data->n, z, data->nabla_F);

  Jac_F_FB(0, data->n, z, F, workV1, workV2, data->nabla_F, H);
}

void FB_compute_error_ncp(void* data_opaque, double* z, double* w, double* Jac_F_merit,
                          double tol, double* err) {
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem*)data_opaque;
  *err = cblas_dnrm2(data->n, Jac_F_merit, 1);
}

void ncp_newton_FBLSA(NonlinearComplementarityProblem* problem, double* z, double* F,
                      int* info, SolverOptions* options) {
  functions_LSA functions_FBLSA_ncp;
  init_lsa_functions(&functions_FBLSA_ncp, &FB_compute_F_ncp, &ncp_FB);
  functions_FBLSA_ncp.compute_H = &FB_compute_H_ncp;
  functions_FBLSA_ncp.compute_error = &FB_compute_error_ncp;

  set_lsa_params_data(options, problem->nabla_F);
  newton_LSA(problem->n, z, F, info, (void*)problem, options, &functions_FBLSA_ncp);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers NCP_NEWTON_FB_FBLSA in the global solver registry.
 */

static void ncp_newton_fblsa_set_default(SolverOptions* options) {
  SOLVER_MAX_ITER(options) = 1000;
  SOLVER_TOL(options) = 1e-4;
}

static int ncp_newton_fblsa_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int ncp_newton_fblsa_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  NonlinearComplementarityProblem* ncp = (NonlinearComplementarityProblem*)problem;
  ncp_newton_FBLSA(ncp, z, w, &info, options);
  return info;
}

static void ncp_newton_fblsa_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_NCP_NEWTON_FB_FBLSA, "NCP_NEWTON_FB_FBLSA",
                       "Newton FBLSA solver for Nonlinear Complementarity Problems",
                       ncp_newton_fblsa_init_wrap, ncp_newton_fblsa_solve_wrap, ncp_newton_fblsa_free_wrap, NULL,
                       ncp_newton_fblsa_set_default, 1000, 1e-4, 0);
