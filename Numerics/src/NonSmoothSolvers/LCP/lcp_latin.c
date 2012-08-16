/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "LA.h"
#include "LCP_Solvers.h"

void lcp_latin(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  /* matrix M of the lcp */
  double * M = problem->M->matrix0;

  /* size of the LCP */
  int n = problem->size;
  int n2 = n * n;


  int i, j,  iter1, nrhs;
  int info2 = 0;
  int itt, it_end;
  int incx, incy;
  int itermax = options->iparam[0];
  double tol = options->dparam[0];
  double k_latin = options->dparam[2];
  double alpha, beta;
  double  err1, num11, err0;
  double res, errmax;
  double den11, den22;
  double  *wc, *zc, *kinvden1, *kinvden2, *wt;
  double *maxwt, *wnum1, *znum1, *ww, *zz;
  double *num1, *kinvnum1, *den1, *den2, *wden1, *zden1;
  double  *kinvwden1, *kzden1;
  double  *k, *kinv, *DPO;

  // char trans='T', notrans='N', uplo='U', diag='N';
  incx = 1;
  incy = 1;

  /* Recup input */


  errmax = tol;
  itt = itermax;


  /* Initialize output */

  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  /* Allocations */

  ww        = (double*) malloc(n * sizeof(double));
  zz        = (double*) malloc(n * sizeof(double));
  wc        = (double*) malloc(n * sizeof(double));
  zc        = (double*) malloc(n * sizeof(double));
  znum1     = (double*) malloc(n * sizeof(double));
  wnum1     = (double*) malloc(n * sizeof(double));
  kinvden1  = (double*) malloc(n * sizeof(double));
  kinvden2  = (double*) malloc(n * sizeof(double));
  wt        = (double*) malloc(n * sizeof(double));
  maxwt     = (double*) malloc(n * sizeof(double));
  num1      = (double*) malloc(n * sizeof(double));
  kinvnum1  = (double*) malloc(n * sizeof(double));
  den1      = (double*) malloc(n * sizeof(double));
  den2      = (double*) malloc(n * sizeof(double));
  wden1     = (double*) malloc(n * sizeof(double));
  zden1     = (double*) malloc(n * sizeof(double));
  kinvwden1 = (double*) malloc(n * sizeof(double));
  kzden1    = (double*) malloc(n * sizeof(double));
  DPO       = (double*) malloc(n2 * sizeof(double));
  k         = (double*) malloc(n2 * sizeof(double));
  kinv      = (double*) malloc(n2 * sizeof(double));



  /* Initialization */

  for (i = 0; i < n2; i++)
  {


    if (i < n)
    {

      wc[i]       = 0.0;
      zc[i]       = 0.0;
      z[i]        = 0.0;
      w[i]        = 0.0;
      znum1[i]    = 0.0;
      wnum1[i]    = 0.0;
      kinvden1[i] = 0.0;
      kinvden2[i] = 0.0;
      wt[i]       = 0.0;
      maxwt[i]    = 0.0;
      num1[i]     = 0.0;
      kinvnum1[i] = 0.0;
      den1[i]     = 0.0;
      den2[i]     = 0.0;
    }

    k[i]          = 0.0;
    kinv[i]       = 0.0;
    DPO[i]        = 0.0;


  }





  for (i = 0 ; i < n ; i++)
  {

    k[i * n + i] =  k_latin * M[i * n + i];

    if (fabs(k[i * n + i]) < DBL_EPSILON)
    {

      if (verbose > 0)
      {
        printf(" Warning nul diagonal term in k matrix \n");
      }

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

      *info = 3;

      return;

    }
    else

      kinv[i + n * i] = 1.0 / k[i + n * i];

  }



  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      DPO[i + n * j] = M[j * n + i] + k[i + n * j];







  /*            Cholesky              */


  DPOTRF(LA_UP, n, DPO , n, &info2);


  if (info2 != 0)
  {
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

    *info = 2;
    return;
  }

  /*            End of Cholesky          */




  /*            Iteration loops  */

  iter1 = 0;
  err1 = 1.;


  while ((iter1 < itt) && (err1 > errmax))
  {

    /*       Linear stage (zc,wc) -> (z,w)*/


    alpha = 1.;
    beta  = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, zc, incx, beta, wc, incy);


    DCOPY(n, problem->q, incx, znum1, incy);


    alpha = -1.;
    DSCAL(n , alpha , znum1 , incx);

    alpha = 1.;
    DAXPY(n, alpha, wc, incx, znum1, incy);
    nrhs = 1;

#ifdef USE_MKL
    DTRTRS(LA_UP, CLA_TRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, info2);
    DTRTRS(LA_UP, CLA_NOTRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, info2);
#else
    DTRTRS(LA_UP, LA_TRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info2);
    DTRTRS(LA_UP, LA_NOTRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info2);
#endif
    DCOPY(n, znum1, incx, z, incy);



    alpha = -1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, z, incx, beta, wc, incy);

    DCOPY(n, wc, incx, w, incy);


    /*         Local Stage                  */

    DCOPY(n, w, incx, wt, incy);
    alpha = -1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, z, incx, beta, wt, incy);

    for (i = 0; i < n; i++)
    {
      if (wt[i] > 0.0)
      {
        wc[i] = wt[i];
        zc[i] = 0.0;
      }
      else
      {
        wc[i] = 0.0;
        zc[i] =  -wt[i] / k[i + n * i];
      }
    }



    /*        Convergence criterium                */


    DCOPY(n, w, incx, wnum1, incy);
    alpha = -1.;
    DAXPY(n, alpha, wc, incx, wnum1, incy);


    DCOPY(n, z, incx, znum1, incy);
    DAXPY(n, alpha, zc, incx, znum1, incy);



    alpha = 1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, znum1, incx, beta, wnum1, incy);


    /*   wnum1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))   */



    alpha = 1.;
    beta = 0.;
    DGEMV(LA_TRANS, n, n, alpha, kinv, n, wnum1, incx, beta, kinvnum1, incy);


    num11 = DDOT(n, wnum1, incx, kinvnum1, incy);




    DCOPY(n, z, incx, zz, incy);
    DCOPY(n, w, incx, ww, incy);

    alpha = 1.;
    DAXPY(n, alpha, wc, incx, ww, incy);

    DAXPY(n, alpha, zc, incx, zz, incy);

    beta = 0.;
    alpha = 1.;
    DGEMV(LA_TRANS, n, n, alpha, k, n, zz, incx, beta, kzden1, incy);


    den22 = DDOT(n, zz, incx, kzden1, incy);

    beta = 0.;
    alpha = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kinv, n, ww, incx, beta, kinvwden1, incy);

    den11 = DDOT(n, ww, incx, kinvwden1, incy);





    err0   = num11 / (den11 + den22);
    err1   = sqrt(err0);

    it_end = iter1;
    res    = err1;

    iter1  = iter1 + 1;

    options->iparam[1] = it_end;
    options->dparam[1] = res;

  }



  if (err1 > errmax)
  {
    if (verbose > 0) printf("No convergence of LATIN after %d iterations, the residue is %g\n", iter1, err1);
    *info = 1;
  }
  else
  {
    if (verbose > 0) printf("Convergence of LATIN after %d iterations, the residue is %g \n", iter1, err1);
    *info = 0;
  }




  free(wc);

  free(DPO);
  free(k);
  free(kinv);

  free(zz);
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

int linearComplementarity_latin_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the Latin Solver\n");
  }


  options->solverId = SICONOS_LCP_LATIN;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->dparam[2] = 0.3;
  options->dparam[3] = 1.0;


  return 0;
}
