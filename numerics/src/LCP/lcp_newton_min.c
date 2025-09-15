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
#include <stdio.h>   // for printf
#include <stdlib.h>  // for free, malloc

#include "LCP_Solvers.h"                   // for lcp_compute_error, lcp_new...
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix
#include "SiconosBlas.h"                   // for cblas_daxpy, cblas_dcopy
#include "SiconosLapack.h"                 // for lapack_int, DGESV
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "numerics_verbose.h"              // for verbose

void lcp_newton_min(LinearComplementarityProblem *problem, double *z, double *w, int *info,
                    SolverOptions *options) {
  /* matrix M/vector q of the lcp */
  assert(problem);
  assert(problem->M);
  double *M = problem->M->matrix0;
  double *q = problem->q;

  /* size of the LCP */
  int n = problem->size;
  assert(n > 0);

  int i, j, iter;
  int m, mm, k;

  int incx, incy;
  double err, a1, b1;
  double alpha;
  lapack_int infoDGESV = 0;

  lapack_int *ipiv;

  double *JacH, *H, *A;

  double *rho;
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tol = options->dparam[SICONOS_DPARAM_TOL];

  incx = 1;
  incy = 1;
  /*input*/

  /*output*/

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  for (i = 0; i < n; i++) {
    z[i] = 1.0;
    w[i] = 1.0;
  }

  /* rho*/
  rho = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; i++) rho[i] = 1.0 / M[i * n + i];
  /* /for (i=0;i<n;i++) rho[i]=1.0/n ;
  // Sizw of the problem*/
  m = 2 * n;
  mm = m * m;
  /* / Creation of the gradient of the function H*/

  JacH = (double *)malloc(m * m * sizeof(double));
  A = (double *)malloc(m * m * sizeof(double));

  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++)
      JacH[j * m + i] = -M[j * n + i]; /* / should be replaced by a tricky use of BLAS*/
  }
  for (j = n; j < m; j++) {
    for (i = 0; i < n; i++) JacH[j * m + i] = 0.0;
    JacH[j * m + j - n] = 1.0;
  }
  for (j = 0; j < m; j++) {
    for (i = n; i < m; i++) JacH[j * m + i] = 0.0;
  }

  /* / Creation of the RHS H, */
  H = (double *)malloc(m * sizeof(double));
  /* / Construction of the RHS*/
  a1 = -1.;
  b1 = -1.;
  /* / q --> H*/
  cblas_dcopy(n, q, incx, H, incy);
  /* / -Mz-q --> H*/
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, a1, M, n, z, incx, b1, H, incy);
  /* / w+H --> H*/
  alpha = 1.0;
  cblas_daxpy(n, alpha, w, incx, H, incy); /* / c'est faux*/

  for (int ii = 0; ii < m - n; ++ii) {
    if (w[ii] > rho[ii] * z[ii]) {
      H[ii + n] = rho[ii] * z[ii];
    } else
      H[ii + n] = w[ii];
  }

  ipiv = (lapack_int *)malloc(m * sizeof(lapack_int));

  iter = 0;
  err = 1.;

  while ((iter < itermax) && (err > tol)) {
    ++iter;
    /* / Construction of the directional derivatives of H, JacH*/
    for (i = 0; i < n; i++) {
      if (w[i] > rho[i] * z[i]) {
        JacH[i * m + i + n] = rho[i];
        JacH[(i + n) * m + (i + n)] = 0.0;
      } else {
        JacH[i * m + i + n] = 0.0;
        JacH[(i + n) * m + (i + n)] = 1.0;
      }
    }

    /* / Computation of the element of the subgradient.*/

    cblas_dcopy(mm, JacH, incx, A, incy);
    k = 1;
    DGESV(m, k, A, m, ipiv, H, m, &infoDGESV);

    if (infoDGESV) {
      if (verbose > 0) {
        printf("Problem in DGESV\n");
      }
      options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
      options->dparam[SICONOS_DPARAM_RESIDU] = err;

      free(H);
      free(A);
      free(JacH);
      free(ipiv);
      free(rho);
      *info = 2;
      return;
    }

    /* / iteration*/
    alpha = -1.0;
    cblas_daxpy(n, alpha, H, incx, z, incy);     /* /  z-H --> z*/
    cblas_daxpy(n, alpha, &H[n], incx, w, incy); /* /  w-H --> w*/

    /* / Construction of the RHS for the next iterate and for the error evalutaion*/
    a1 = 1.;
    b1 = 1.;
    cblas_dcopy(n, q, incx, H, incy); /* / q --> H*/
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, a1, M, n, z, incx, b1, H,
                incy); /* / Mz+q --> H*/
    alpha = -1.0;
    cblas_daxpy(n, alpha, w, incx, H, incy); /* / w-Mz-q --> H*/

    for (i = 0; i < m - n; i++) {
      if (w[i] > rho[i] * z[i])
        H[i + n] = rho[i] * z[i];
      else
        H[i + n] = w[i];
    }

    /* **** Criterium convergence **** */
    lcp_compute_error(problem, z, w, tol, &err);
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if (err > tol) {
    printf(
        "Siconos/Numerics: lcp_newton_min: No convergence of NEWTON_MIN after %d iterations\n",
        iter);
    printf("Siconos/Numerics: lcp_newton_min: The residue is : %g \n", err);
    *info = 1;
  } else {
    if (verbose > 0) {
      printf(
          "Siconos/Numerics: lcp_newton_min: Convergence of NEWTON_MIN after %d iterations\n",
          iter);
      printf("Siconos/Numerics: lcp_newton_min: The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(H);
  free(A);
  free(JacH);
  free(ipiv);
  free(rho);
}
