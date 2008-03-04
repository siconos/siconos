/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include "Relay_Solvers.h"

void pr_latin(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{

  double* vec = problem->M->matrix0;
  double* qq = problem->q;
  int n = problem -> size;
  double *a = problem->a;
  double *b = problem->b;

  double k_latin = options->dparam[2];
  int itt = options->iparam[0];
  double errmax = options->dparam[0];

  int i, j, iter1, nrhs, info2;

  int incx = 1, incy = 1;
  double alpha, beta, mina;
  double  err1, num11, err0;
  double den11, den22;
  double *wc, *zc, *wnum1, *znum1;
  double *zt;
  double *k, *kinv, *DPO;

  /*  char trans='T',uplo='U',notrans='N',diag='N'; */




  /*             Allocations                           */

  k     = (double*)malloc(n * n * sizeof(double));
  DPO   = (double*)malloc(n * n * sizeof(double));
  kinv  = (double*)malloc(n * n * sizeof(double));
  wc    = (double*) malloc(n * sizeof(double));
  zc    = (double*) malloc(n * sizeof(double));
  znum1 = (double*) malloc(n * sizeof(double));
  wnum1 = (double*) malloc(n * sizeof(double));
  zt    = (double*) malloc(n * sizeof(double));


  /*             Initialization                          */


  for (i = 0; i < n * n; i++)
  {
    k[i]    = 0.0;
    kinv[i] = 0.0;
    DPO[i]  = 0.0;
    if (i < n)
    {
      wc[i]    = 0.0;
      zc[i]    = 0.0;
      z[i]     = 0.0;
      w[i]     = 0.0;
      znum1[i] = 0.0;
      wnum1[i] = 0.0;
      zt[i]    = 0.0;
    }
  }


  for (i = 0; i < n; i++)
  {
    k[i + n * i] =  k_latin * vec[i * n + i];

    if (fabs(k[i + n * i]) < 1e-12)
    {

      if (verbose > 0)
        printf("\n Warning nul diagonal term in k matrix \n ");

      free(k);
      free(kinv);
      free(DPO);
      free(wc);
      free(zc);
      free(znum1);
      free(wnum1);
      free(zt);

      *info = 3;
      return;

    }
    else

      kinv[i + n * i] = 1 / k[i + n * i];
  }






  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      DPO[i + n * j] =  vec[j * n + i] + k[i + n * j];
    }


  /*             Cholesky                           */

  DPOTRF(LA_UP, n, DPO, n, info2);



  if (info2 != 0)
  {
    printf("\n Matter with Cholesky factorization \n");

    free(k);
    free(kinv);
    free(DPO);
    free(wc);
    free(zc);
    free(znum1);
    free(wnum1);
    free(zt);

    *info = 2;
    return;
  }


  /*             End of Cholesky                        */



  /*              Iteration loops                         */

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {

    /*          Linear stage (zc,wc) -> (z,w)            */

    alpha = 1.;
    beta  = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, zc, incx, beta, wc, incy);

    DCOPY(n, qq, incx, znum1, incy);


    alpha = -1.;
    DSCAL(n , alpha , znum1 , incx);

    alpha = 1.;
    DAXPY(n, alpha, wc, incx, znum1, incy);


    nrhs = 1;
    DTRTRS(LA_UP, LA_TRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, info2);

    DTRTRS(LA_UP, LA_NOTRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, info2);

    DCOPY(n, znum1, incx, z, incy);

    DCOPY(n, wc, incx, w, incy);

    alpha = -1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, z, incx, beta, w, incy);

    /*            Local stage (z,w)->(zc,wc)            */


    DCOPY(n, z, incx, zt, incy);

    alpha = -1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kinv, n, w, incx, beta, zt, incy);

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
      if (-b[i] < mina)
      {
        zc[i] = mina;
      }
      else
      {
        zc[i] = -b[i];
      }
    }

    DCOPY(n, w, incx, wc, incy);

    alpha = -1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, z, incx, beta, wc, incy);

    alpha = 1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, zc, incx, beta, wc, incy);

    /*             Convergence criterium                     */

    DCOPY(n, w, incx, wnum1, incy);
    alpha = -1.;
    DAXPY(n, alpha, wc, incx, wnum1, incy);
    DCOPY(n, z, incx, znum1, incy);
    DAXPY(n, alpha, zc, incx, znum1, incy);
    alpha = 1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, znum1, incx, beta, wnum1, incy);

    /*    num1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))   */

    num11 = 0.;
    alpha = 1.;
    beta = 0.;
    DGEMV(LA_TRANS, n, n, alpha, kinv, n, wnum1, incx, beta, znum1, incy);

    num11 = DDOT(n, wnum1, incx, znum1, incy);

    DCOPY(n, z, incx, znum1, incy);

    DAXPY(n, alpha, zc, incx, znum1, incy);

    beta = 0.;
    alpha = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, znum1, incx, beta, wnum1, incy);

    den22 = DDOT(n, znum1, incx, wnum1, incy);

    DCOPY(n, w, incx, wnum1, incy);

    alpha = 1.;
    DAXPY(n, alpha, wc, incx, wnum1, incy);

    beta = 0.;
    alpha = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kinv, n, wnum1, incx, beta, znum1, incy);

    den11 = DDOT(n, wnum1, incx, znum1, incy);

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
  free(zt);
  free(DPO);
  free(k);
  free(kinv);



}
