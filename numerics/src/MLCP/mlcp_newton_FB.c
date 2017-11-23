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
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include "SiconosLapack.h"
#include "Newton_methods.h"
#include "FischerBurmeister.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

static void FB_compute_F_mlcp(void* data_opaque, double* z, double* w)
{
  // Computation of the new value w = F(z) = Mz + q
  // q --> w
  MixedLinearComplementarityProblem* problem = (MixedLinearComplementarityProblem *)data_opaque;
  assert(problem->M);
  assert(problem->M->matrix0);
  unsigned int n = problem->n;
  unsigned int m = problem->m;
  /* Problem in the form (M,q) */
  if (problem->isStorageType1)
  {
    cblas_dcopy(n, problem->q, 1, w, 1);
    // Mz+q --> w
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, problem->M->matrix0, n, z, 1, 1.0, w, 1);
  }
  else
  {
    /* Links to problem data */
    double *a = problem->q;
    double *b = &problem->q[n];
    double *A = problem->A;
    double *B = problem->B;
    double *C = problem->C;
    double *D = problem->D;

    /* Compute "equalities" part, we = Au + Cv + a - Must be equal to 0 */
    cblas_dcopy(n, a , 1 , w , 1); //  we = w[0..n-1] <-- a
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, A, n, z, 1, 1.0, w, 1); // we <-- A*u + we
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1.0, C, m, &z[n], 1, 1.0, w, 1); // we <-- C*v + we

    /* Computes part which corresponds to complementarity */
    double* w_c = &w[n]; // No copy!!
    cblas_dcopy(m, b, 1, w_c, 1); //  wi = w[n..m] <-- b
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, D, n, z, 1, 1.0, w_c, 1); // we <-- D*u + we
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, m, 1.0, B, m, &z[n], 1, 1.0, w_c, 1); // we <-- B*v + we
  }
 }

static void FB_compute_H_mlcp(void* data_opaque, double* z, double* w, double* workV1, double* workV2, NumericsMatrix* H)
{
  printf("MLCP FB_compute_H_mlcp not implemented yet");
  exit(1);
#if 0
  MixedLinearComplementarityProblem* data = (MixedLinearComplementarityProblem *)data_opaque;
  unsigned int n = data->size;
  assert(data->M);
  assert(data->M->matrix0);
  double* M = data->M->matrix0;
  double normi;

  // workV1 = "z" in Facchibei--Pang p. 808
  // "z_i" = 1 if z_i = w_i = 0.0
  // M^T.workV1 --> workV2
  cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, M, n , workV1, 1, 0.0, workV2, 1);
  for (unsigned int i = 0; i < n; ++i)
  {
    if (workV1[i] != 0.0) // i in beta
    {
      normi = sqrt(workV1[i] * workV1[i] + workV2[i] * workV2[i]);
      for (unsigned int j = 0; j < n; j++)
      {
        H[j * n + i] = (workV2[i] / normi - 1.0) * M[j * n + i];
      }
      H[i * n + i] += (workV1[i] / normi - 1.0);

    }
    else // i not in beta
    {
      normi = sqrt(z[i] * z[i] + w[i] * w[i]);
      for (unsigned int j = 0; j < n; j++)
      {
        H[j * n + i] = (w[i] / normi - 1.0) * M[j * n + i];
      }
      H[i * n + i] += (z[i] / normi - 1.0);
    }

  }
#endif
}

void mlcp_mixed_FB(void* data_opaque, double* z, double* F, double* F_FB)
{
  phi_Mixed_FB(((MixedLinearComplementarityProblem *)data_opaque)->n, ((MixedLinearComplementarityProblem *)data_opaque)->m, z, F, F_FB);
}

void FB_compute_error_mlcp(void* data_opaque, double* z, double* w, double* nabla_theta, double tol, double* err)
{
  mlcp_compute_error((MixedLinearComplementarityProblem *)data_opaque, z, w, tol, err);
}

void mlcp_newton_FB(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  functions_LSA functions_FBLSA_mlcp;
  init_lsa_functions(&functions_FBLSA_mlcp, &FB_compute_F_mlcp, (compute_F_merit_ptr)&mlcp_FB);
  functions_FBLSA_mlcp.compute_H = &FB_compute_H_mlcp;
  functions_FBLSA_mlcp.compute_error = &FB_compute_error_mlcp;

  set_lsa_params_data(options, problem->M);
  newton_LSA(problem->n + problem->m, z, w, info, (void *)problem, options, &functions_FBLSA_mlcp);
}


