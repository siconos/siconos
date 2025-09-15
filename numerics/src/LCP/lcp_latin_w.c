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

#include <math.h>    // for fabs
#include <stdio.h>   // for printf, NULL
#include <stdlib.h>  // for free, malloc

#include "LCP_Solvers.h"                   // for lcp_compute_error_only
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_daxpy
#include "SiconosLapack.h"     // for DTRTRS, DPOTRF, lapack_int, LA_UP, LA_NONUNIT, LA_NOTRANS
#include "SolverOptions.h"     // for SolverOptions, solver_opti...
#include "lcp_cst.h"           // for SICONOS_LCP_DPARAM_LATIN_P...
#include "numerics_verbose.h"  // for verbose

void lcp_latin_w(LinearComplementarityProblem *problem, double *z, double *w, int *info,
                 SolverOptions *options) {
  /* matrix M/vector q of the lcp */
  double *M = problem->M->matrix0;

  double *q = problem->q;

  /* size of the LCP */
  int n = problem->size;
  int n2 = n * n;

  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  double k_latin = options->dparam[SICONOS_LCP_DPARAM_LATIN_PARAMETER];
  double omega = options->dparam[SICONOS_LCP_DPARAM_RHO];

  int i, j, iter1, nrhs;
  lapack_int info2 = 0;
  int itt, it_end;
  int incx, incy;

  double alpha, beta;
  double err1;
  double res, errmax;
  double *wc, *zc, *kinvden1, *kinvden2, *wt;
  double *maxwt, *wnum1, *znum1, *ww, *zz;
  double *num1, *kinvnum1, *den1, *den2, *wden1, *zden1;
  double *kinvwden1, *kzden1;
  double *k, *kinv, *DPO;
  double *zn, *wn;

  /* char trans='T', notrans='N', uplo='U', diag='N'; */

  incx = 1;
  incy = 1;

  errmax = tol;
  itt = itermax;

  /* Initialize output */

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  /* Allocations */

  ww = (double *)malloc(n * sizeof(double));
  zz = (double *)malloc(n * sizeof(double));
  wc = (double *)malloc(n * sizeof(double));
  zc = (double *)malloc(n * sizeof(double));
  znum1 = (double *)malloc(n * sizeof(double));
  wnum1 = (double *)malloc(n * sizeof(double));
  kinvden1 = (double *)malloc(n * sizeof(double));
  kinvden2 = (double *)malloc(n * sizeof(double));
  wt = (double *)malloc(n * sizeof(double));
  maxwt = (double *)malloc(n * sizeof(double));
  num1 = (double *)malloc(n * sizeof(double));
  kinvnum1 = (double *)malloc(n * sizeof(double));
  den1 = (double *)malloc(n * sizeof(double));
  den2 = (double *)malloc(n * sizeof(double));
  wden1 = (double *)malloc(n * sizeof(double));
  zden1 = (double *)malloc(n * sizeof(double));
  kinvwden1 = (double *)malloc(n * sizeof(double));
  kzden1 = (double *)malloc(n * sizeof(double));
  zn = (double *)malloc(n * sizeof(double));
  wn = (double *)malloc(n * sizeof(double));

  DPO = (double *)malloc(n2 * sizeof(double));
  k = (double *)malloc(n2 * sizeof(double));
  kinv = (double *)malloc(n2 * sizeof(double));

  /* Initialization */

  for (int i = 0; i < n2; i++) {
    if (i < n) {
      wc[i] = 0.0;
      zc[i] = 0.0;
      z[i] = 0.0;
      zn[i] = 0.0;
      wn[i] = 0.0;
      w[i] = 0.0;
      znum1[i] = 0.0;
      wnum1[i] = 0.0;
      kinvden1[i] = 0.0;
      kinvden2[i] = 0.0;
      wt[i] = 0.0;
      maxwt[i] = 0.0;
      num1[i] = 0.0;
      kinvnum1[i] = 0.0;
      den1[i] = 0.0;
      den2[i] = 0.0;
    }

    k[i] = 0.0;
    kinv[i] = 0.0;
    DPO[i] = 0.0;
  }

  for (i = 0; i < n; i++) {
    k[i * n + i] = k_latin * M[i * n + i];

    if (fabs(k[i * n + i]) < 1e-12) {
      if (verbose > 0) {
        printf(" Warning nul diagonal term in k matrix \n");
      }

      free(ww);

      free(zz);
      free(wn);
      free(zn);
      free(wc);
      free(zc);
      free(znum1);
      free(wnum1);
      free(kinvden1);
      free(kinvden2);
      free(wt);
      free(maxwt);
      free(num1);
      free(kinvnum1);
      free(den1);
      free(den2);
      free(wden1);
      free(zden1);
      free(kinvwden1);
      free(kzden1);
      free(DPO);
      free(k);
      free(kinv);

      *info = 3;

      return;

    } else

      kinv[i + n * i] = 1.0 / k[i + n * i];
  }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) DPO[i + n * j] = M[j * n + i] + k[i + n * j];

  /*            Cholesky              */

  DPOTRF(LA_UP, n, DPO, n, &info2);

  if (info2 != 0) {
    printf(" Matter with Cholesky Factorization \n ");

    free(ww);
    free(zz);
    free(wc);
    free(zc);
    free(znum1);
    free(wnum1);
    free(kinvden1);
    free(kinvden2);
    free(wt);
    free(maxwt);
    free(num1);
    free(kinvnum1);
    free(den1);
    free(den2);
    free(wden1);
    free(zden1);
    free(kinvwden1);
    free(kzden1);
    free(DPO);
    free(k);
    free(kinv);
    free(zn);
    free(wn);

    *info = 2;
    return;
  }

  /*            End of Cholesky          */

  /*            Iteration loops  */

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax)) {
    /*       Linear stage (zc,wc) -> (z,w)*/

    alpha = 1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, k, n, zc, incx, beta, wc, incy);

    cblas_dcopy(n, q, incx, znum1, incy);

    alpha = -1.;
    cblas_dscal(n, alpha, znum1, incx);

    alpha = 1.;
    cblas_daxpy(n, alpha, wc, incx, znum1, incy);
    nrhs = 1;

    DTRTRS(LA_UP, LA_TRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info2);
    DTRTRS(LA_UP, LA_NOTRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info2);
    cblas_dcopy(n, znum1, incx, z, incy);

    alpha = -1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, k, n, z, incx, beta, wc, incy);

    cblas_dcopy(n, wc, incx, w, incy);

    alpha = omega;
    cblas_dscal(n, alpha, z, incx);

    alpha = 1.0 - omega;
    cblas_daxpy(n, alpha, zn, incx, z, incy);

    alpha = omega;
    cblas_dscal(n, alpha, w, incx);

    alpha = 1.0 - omega;
    cblas_daxpy(n, alpha, wn, incx, w, incy);

    cblas_dcopy(n, w, incx, wn, incy);
    cblas_dcopy(n, z, incx, zn, incy);

    /*         Local Stage                  */

    cblas_dcopy(n, w, incx, wt, incy);
    alpha = -1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, k, n, z, incx, beta, wt, incy);

    for (i = 0; i < n; i++) {
      if (wt[i] > 0.0) {
        wc[i] = wt[i];
        zc[i] = 0.0;
      } else {
        wc[i] = 0.0;
        zc[i] = -wt[i] / k[i + n * i];
      }
    }

    /*        Convergence criterium                */

    cblas_dcopy(n, w, incx, wnum1, incy);
    alpha = -1.;
    cblas_daxpy(n, alpha, wc, incx, wnum1, incy);

    cblas_dcopy(n, z, incx, znum1, incy);
    cblas_daxpy(n, alpha, zc, incx, znum1, incy);

    alpha = 1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, k, n, znum1, incx, beta, wnum1, incy);

    /*   wnum1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))   */

    alpha = 1.;
    beta = 0.;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, kinv, n, wnum1, incx, beta, kinvnum1,
                incy);

    cblas_dcopy(n, z, incx, zz, incy);
    cblas_dcopy(n, w, incx, ww, incy);

    alpha = 1.;
    cblas_daxpy(n, alpha, wc, incx, ww, incy);

    cblas_daxpy(n, alpha, zc, incx, zz, incy);

    beta = 0.;
    alpha = 1.;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, k, n, zz, incx, beta, kzden1, incy);

    beta = 0.;
    alpha = 1.;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, alpha, kinv, n, ww, incx, beta, kinvwden1,
                incy);

    lcp_compute_error_only(n, z, w, &err1);

    it_end = iter1;
    res = err1;

    iter1 = iter1 + 1;

    options->iparam[SICONOS_IPARAM_ITER_DONE] = it_end;
    options->dparam[SICONOS_DPARAM_RESIDU] = res;
  }

  if (err1 > errmax) {
    if (verbose > 0)
      printf("No convergence of LATIN_W after %d iterations, the residue is %g\n", iter1,
             err1);
    *info = 1;
  } else {
    if (verbose > 0)
      printf("Convergence of LATIN_W after %d iterations, the residue is %g \n", iter1, err1);
    *info = 0;
  }

  free(wc);

  free(DPO);
  free(k);
  free(kinv);

  free(zz);
  free(zn);
  free(wn);
  free(ww);
  free(zc);
  free(znum1);
  free(wnum1);
  free(kinvden1);
  free(kinvden2);
  free(wt);
  free(maxwt);
  free(num1);
  free(kinvnum1);
  free(den1);
  free(den2);
  free(wden1);
  free(zden1);
  free(kinvwden1);
  free(kzden1);
}
void lcp_latin_w_set_default(SolverOptions *options) {
  options->dparam[SICONOS_LCP_DPARAM_LATIN_PARAMETER] = 0.3;
  options->dparam[SICONOS_LCP_DPARAM_RHO] = 1.0;
}
