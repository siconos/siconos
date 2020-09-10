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

#include "NCP_Solvers.h"                      // for ncp_newton_minFBLSA
#include "Newton_methods.h"                   // for functions_LSA, init_lsa...
#include "NonlinearComplementarityProblem.h"  // for NonlinearComplementarit...
#include "NumericsFwd.h"                      // for NonlinearComplementarit...
#include "min_merit.h"                        // for F_min, Jac_F_min
#include "ncp_newton_FBLSA.h"                 // for FB_compute_F_ncp, FB_co...

static void ncp_min(void* data_opaque, double* z, double* F, double* Fmin)
{
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem *)data_opaque;

  F_min(0, data->n, z, F, Fmin);
}

static void min_compute_H_ncp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, NumericsMatrix* H)
{
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem *)data_opaque;

  data->compute_nabla_F(data->env, data->n, z, data->nabla_F);

  Jac_F_min(0, data->n, z, F, data->nabla_F, H);
}

void ncp_newton_minFBLSA(NonlinearComplementarityProblem* problem, double *z, double* F, int *info, SolverOptions* options)
{
  functions_LSA functions_minFBLSA_ncp;
  init_lsa_functions(&functions_minFBLSA_ncp, &FB_compute_F_ncp, &ncp_FB);
  functions_minFBLSA_ncp.compute_H = &FB_compute_H_ncp;
  functions_minFBLSA_ncp.compute_error = &FB_compute_error_ncp;
  functions_minFBLSA_ncp.compute_RHS_desc = &ncp_min;
  functions_minFBLSA_ncp.compute_H_desc = &min_compute_H_ncp;

  set_lsa_params_data(options, problem->nabla_F);
  newton_LSA(problem->n, z, F, info, (void *)problem, options, &functions_minFBLSA_ncp);
}
