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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "Relay_Solvers.h"
#include "SiconosLapack.h"

void dr_latin(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{

  double* vec = problem->M->matrix0;
  double* qq = problem->q;
  int n = problem -> size;
  double *a = problem->ub;
  double *b = problem->lb;
  //\todo Rewrite completely the algorithm with a projection.
  int ib;
  for (ib = 0; ib < n; ib++) b[ib] = -b[ib];

  double k_latin = options->dparam[2];
  int itermax = options->iparam[0];
  double errmax = options->dparam[0];

  int i, j, iter1, nrhs;
  int info2 = 0;
  int incx = 1, incy = 1;

  double alpha, beta, mina;
  double err1, num11, err0;
  double den11, den22;
  double *wc, *zc, *wt, *wnum1, *znum1;
  double *zt, *kinvnum1;

  /*  char trans='T',notrans='N', uplo='U', diag='N'; */

  double *k, *kinv, *DPO;




  /*             Allocations                           */

  k         = (double*) malloc(n * n * sizeof(double));
  DPO       = (double*) malloc(n * n * sizeof(double));
  kinv      = (double*) malloc(n * n * sizeof(double));
  wc        = (double*) malloc(n * sizeof(double));
  zc        = (double*) malloc(n * sizeof(double));
  znum1     = (double*) malloc(n * sizeof(double));
  wnum1     = (double*) malloc(n * sizeof(double));
  wt        = (double*) malloc(n * sizeof(double));
  zt        = (double*) malloc(n * sizeof(double));
  kinvnum1  = (double*) malloc(n * sizeof(double));


  /*             Initialisation                   */


  for (i = 0; i < n * n; i++)
  {
    k[i]    =  0.0;
    DPO[i]  =  0.0;
    kinv[i] =  0.0;

    if (i < n)
    {
      wc[i]       = 0.0;
      zc[i]       = 0.0;
      z[i]        = 0.0;
      w[i]        = 0.0;
      znum1[i]    = 0.0;
      wnum1[i]    = 0.0;
      wt[i]       = 0.0;
      zt[i]       = 0.0;
      kinvnum1[i] = 0.0;
    }

  }


  for (i = 0; i < n; i++)
  {
    k[i + n * i] = k_latin * vec[i * n + i];

    if (fabs(k[i + n * i]) < DBL_EPSILON)
    {

      if (verbose > 0)
        printf("\n Warning nul diagonal term in k matrix \n");

      free(k);
      free(kinv);
      free(DPO);
      free(wc);
      free(zc);
      free(znum1);
      free(wnum1);
      free(wt);
      free(zt);
      free(kinvnum1);

      *info = 3;
      return;

    }
    else

      kinv[i + n * i] = 1 / k[i + n * i];

  }



  for (j = 0; j < n; j++)
  {
    for (i = 0; i < n; i++)

    {
      DPO[i + n * j] =  vec[j * n + i] + k[i + n * j];
    }
  }



  /*                    Cholesky             */



  DPOTRF(LA_UP, n, DPO , n, &info2);

  if (info2 != 0)
  {

    if (verbose > 0)
      printf("\n Matter with Cholesky factorization \n");

    free(k);
    free(kinv);
    free(DPO);
    free(wc);
    free(zc);
    free(znum1);
    free(wnum1);
    free(wt);
    free(zt);
    free(kinvnum1);

    *info = 2;
    return;
  }



  /*            End of cholesky                   */






  /*             Iteration loops        */


  iter1  = 0;
  err1   = 1.;


  while ((iter1 < itermax) && (err1 > errmax))
  {

    /*   Linear stage (zc,wc) -> (z,w)       */


    alpha = 1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, k, n, zc, incx, beta, wc, incy);

    cblas_dcopy(n, qq, incx, znum1, incy);


    alpha = -1.;
    cblas_dscal(n, alpha, znum1, incx);

    alpha = 1.;
    cblas_daxpy(n, alpha, wc, incx, znum1, incy);

    nrhs = 1;
    DTRTRS(LA_UP, LA_TRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info2);
    DTRTRS(LA_UP, LA_NOTRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info2);

    cblas_dcopy(n, znum1, incx, z, incy);

    cblas_dcopy(n, wc, incx, w, incy);

    alpha = -1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, k, n, z, incx, beta, w, incy);

    /*     Local stage (z,w)->(zc,wc)         */


    cblas_dcopy(n, w, incx, zt, incy);

    alpha = -1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, k, n, z, incx, beta, zt, incy);

    for (i = 0; i < n; i++)
    {
      if (a[i] < zt[i])
      {
        mina = a[i];
      }
      else
      {
        mina = zt[i];
      }
      if (mina > -b[i])
      {
        wc[i] = mina;
      }
      else
      {
        wc[i] = -b[i];
      }
    }

    cblas_dcopy(n, wc, incx, wnum1, incy);

    alpha = -1.;
    cblas_daxpy(n, alpha, zt, incx, wnum1, incy);

    cblas_dcopy(n, wnum1, incx, zt, incy);

    alpha = 1.;
    beta = 0.;
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, kinv, n, zt, incx, beta, zc, incy);

    /*            Convergence criterium          */

    cblas_dcopy(n, w, incx, wnum1, incy);

    alpha = -1.;
    cblas_daxpy(n, alpha, wc, incx, wnum1, incy);

    cblas_dcopy(n, z, incx, znum1, incy);

    cblas_daxpy(n, alpha, zc, incx, znum1, incy);

    alpha = 1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, k, n, znum1, incx, beta, wnum1, incy);

    /*       num1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))  */

    num11 = 0.;
    alpha = 1.;
    beta = 0.;
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, kinv, n, wnum1, incx, beta, kinvnum1, incy);

    num11 = cblas_ddot(n, wnum1, incx, kinvnum1, incy);

    cblas_dcopy(n, w, incx, wnum1, incy);

    alpha = 1.;
    cblas_daxpy(n, alpha, wc, incx, wnum1, incy);

    cblas_dcopy(n, z, incx, znum1, incy);

    cblas_daxpy(n, alpha, zc, incx, znum1, incy);

    beta = 0.;
    alpha = 1.;
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, k, n, znum1, incx, beta, kinvnum1, incy);

    den22 = cblas_ddot(n, znum1, incx, kinvnum1, incy);

    beta = 0.;
    alpha = 1.;

    cblas_dgemv(CblasColMajor,CblasTrans, n, n, alpha, kinv, n, wnum1, incx, beta, kinvnum1, incy);

    den11 = cblas_ddot(n, wnum1, incx, kinvnum1, incy);

    err0 = num11 / (den11 + den22);
    err1 = sqrt(err0);
    iter1 = iter1 + 1;


    options->iparam[1] = iter1;
    options->dparam[1] = err1;



  }


  if (err1 > errmax)
  {
    if (verbose > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter1, err1);

    *info = 1;
  }
  else
  {
    if (verbose > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter1, err1);
    *info = 0;
  }


  free(wc);
  free(zc);
  free(znum1);
  free(wnum1);
  free(kinvnum1);
  free(wt);
  free(zt);

  free(k);
  free(DPO);
  free(kinv);
}
