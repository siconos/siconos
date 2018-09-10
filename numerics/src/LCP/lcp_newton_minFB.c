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

#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"

#include "SiconosLapack.h"
#include "Newton_methods.h"
#include "FischerBurmeister.h"
#include "min_merit.h"
#include "lcp_newton_FB.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

static void lcp_min(void* data_opaque, double* z, double* F, double* Fmin)
{
  F_min(0, ((LinearComplementarityProblem *)data_opaque)->size, z, F, Fmin);
}

static void min_compute_H_lcp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, NumericsMatrix* H)
{
  LinearComplementarityProblem* data = (LinearComplementarityProblem *)data_opaque;
  unsigned int n = data->size;
  assert(data->M);

  Jac_F_min(0, n, z, F, data->M, H);
}

void lcp_newton_minFB(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  functions_LSA functions_minFBLSA_lcp;
  init_lsa_functions(&functions_minFBLSA_lcp, &FB_compute_F_lcp, &lcp_FB);
  functions_minFBLSA_lcp.compute_H = &FB_compute_H_lcp;
  functions_minFBLSA_lcp.compute_error = &FB_compute_error_lcp;
  functions_minFBLSA_lcp.compute_RHS_desc = &lcp_min;
  functions_minFBLSA_lcp.compute_H_desc = &min_compute_H_lcp;

  set_lsa_params_data(options, problem->M);
  newton_LSA(problem->size, z, w, info, (void *)problem, options, &functions_minFBLSA_lcp);
}
