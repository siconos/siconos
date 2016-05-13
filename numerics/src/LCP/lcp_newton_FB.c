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

#include "LCP_Solvers.h"
#include "SiconosLapack.h"
#include "Newton_Methods.h"
#include "FischerBurmeister.h"

#include "lcp_newton_FB.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void FB_compute_F_lcp(void* data_opaque, double* z, double* w)
{
  // Computation of the new value w = F(z) = Mz + q
  // q --> w
  LinearComplementarityProblem* data = (LinearComplementarityProblem *)data_opaque;
  assert(data->M);
  assert(data->M->matrix0);
  unsigned int n = data->size;
  cblas_dcopy(n , data->q, 1, w, 1);
  // Mz+q --> w
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, data->M->matrix0, n, z, 1, 1.0, w, 1);
}

void FB_compute_H_lcp(void* data_opaque, double* z, double* w, double* workV1, double* workV2, NumericsMatrix* H)
{
  LinearComplementarityProblem* data = (LinearComplementarityProblem *)data_opaque;
  unsigned int n = data->size;
  assert(data->M);

  Jac_F_FB(0, n, z, w, workV1, workV2, data->M, H);
}

void FB_compute_error_lcp(void* data_opaque, double* z, double* w, double* notused, double tol, double* err)
{
  lcp_compute_error((LinearComplementarityProblem *)data_opaque, z, w, tol, err);
}

void lcp_FB(void* data_opaque, double* z, double* F, double* F_FB)
{
  phi_FB(((LinearComplementarityProblem *)data_opaque)->size, z, F, F_FB);
}

void lcp_newton_FB(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  functions_LSA functions_FBLSA_lcp;
  init_lsa_functions(&functions_FBLSA_lcp, &FB_compute_F_lcp, &lcp_FB);
  functions_FBLSA_lcp.compute_H = &FB_compute_H_lcp;
  functions_FBLSA_lcp.compute_error = &FB_compute_error_lcp;

  set_lsa_params_data(options, problem->M);
  newton_LSA(problem->size, z, w, info, (void *)problem, options, &functions_FBLSA_lcp);
}
