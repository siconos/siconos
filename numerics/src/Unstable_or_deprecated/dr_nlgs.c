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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Relay_Solvers.h"
#include "SiconosBlas.h"

void dr_nlgs(RelayProblem *problem, double *z, double *w, int *info, SolverOptions *options) {
  double *vec = problem->M->matrix0;
  double *q = problem->q;
  int n = problem->size;
  double *a = problem->ub;
  double *b = problem->lb;
  //\todo Rewrite completely the algorithm with a projection.
  int ib;
  for (ib = 0; ib < n; ib++) b[ib] = -b[ib];
  int itt = options->iparam[0];
  double errmax = options->dparam[0];

  int i, j, iter1, k;
  int incx = 1, incy = 1;
  double alpha, beta, mina;
  double err1, num, den, avn, xn, apn;
  double *zt, *wnum1;

  wnum1 = (double *)malloc(n * sizeof(double));
  zt = (double *)malloc(n * sizeof(double));

  for (i = 0; i < n; i++) {
    w[i] = 0.;
    z[i] = 0.;
    zt[i] = 0.;
    wnum1[i] = 0.;
  }

  iter1 = 1;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax)) {
    iter1 = iter1 + 1;

    for (i = 0; i < n; i++) {
      avn = 0.;
      apn = 0.;

      for (j = 0; j <= i - 1; j++) avn = avn + vec[j * n + i] * z[j];

      for (k = i + 1; k < n; k++) apn = apn + vec[k * n + i] * z[k];

      xn = -q[i] - avn - apn;

      zt[i] = -xn;

      if (a[i] < zt[i]) {
        mina = a[i];
      } else {
        mina = zt[i];
      }

      if (-b[i] < mina) {
        w[i] = mina;
      } else {
        w[i] = -b[i];
      }

      if (fabs(vec[i * n + i]) < 1e-12) {
        printf("\n Warning nul diagonal term of M \n");

        free(zt);
        free(wnum1);

        *info = 2;

        return;

      } else
        z[i] = 1 / vec[i * n + i] * (w[i] + xn);
    }

    /*              Convergence criterium              */

    cblas_dcopy(n, w, incx, wnum1, incy);

    alpha = -1.;
    cblas_daxpy(n, alpha, q, incx, wnum1, incy);

    alpha = 1.;
    beta = -1.;
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, vec, n, z, incx, beta, wnum1, incy);

    num = cblas_ddot(n, wnum1, incx, wnum1, incy);

    den = cblas_ddot(n, q, incx, q, incy);

    err1 = sqrt(num) / sqrt(den);
    options->iparam[1] = iter1;
    options->dparam[1] = err1;
  }

  if (err1 > errmax) {
    if (verbose > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter1, err1);

    *info = 1;
  } else {
    if (verbose > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter1, err1);

    *info = 0;
  }

  free(wnum1);
  free(zt);
}
