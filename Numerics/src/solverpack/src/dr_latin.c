/*!\file dr_latin.c
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"


void dr_latin(double *vec , double *qq , int *nn , double *k_latin , double *a , double *b , int *itermax , double * tol ,
              double *z , double *w , int *it_end , double *res , int *info)
{

  FILE *f101;
  int i, j, kk, iter1, zzz;
  int n = *nn, incx = 1, incy = 1;
  double errmax = *tol, alpha, beta, mina, maxa, bb, cc, zw, aa;
  double rr, rrr, r1, r2, r3, invR0, invRT0, err1, z0, num11, err0, invRTinvR0;
  double den11, den22, vv;
  double *wc, *zc, *wt, *wnum1, *znum1;
  double *zt, *kinvnum1;
  char trans = 'T';
  double **k, **A, **R, **RT, **invRT, **invR;
  double **invRTinvR, **kinv, xpivot;

  f101 = fopen("resultat_latin.dat", "w+");

  A         = (double **)malloc(n * sizeof(double*));
  k         = (double **)malloc(n * sizeof(double*));
  R         = (double **)malloc(n * sizeof(double*));
  RT        = (double **)malloc(n * sizeof(double*));
  invR      = (double **)malloc(n * sizeof(double*));
  kinv      = (double **)malloc(n * sizeof(double*));
  invRT     = (double **)malloc(n * sizeof(double*));
  invRTinvR = (double **)malloc(n * sizeof(double*));

  for (i = 0 ; i < n ; ++i)
  {
    A[i]      = (double *)malloc(n * sizeof(double));
    k[i]      = (double *)malloc(n * sizeof(double));
    R[i]      = (double *)malloc(n * sizeof(double));
    RT[i]     = (double *)malloc(n * sizeof(double));
    invR[i]   = (double *)malloc(n * sizeof(double));
    kinv[i]   = (double *)malloc(n * sizeof(double));
    invRT[i]  = (double *)malloc(n * sizeof(double));
    invRTinvR = (double *)malloc(n * sizeof(double));
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      k[i][j] = 0.;
    }
  }

  for (i = 0; i < n; i++)
  {
    k[i][i] = *k_latin * vec[i * n + i];
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      A[i][j] =  vec[i * n + j] + k[i][j];
      R[i][j] =  0.;
    }
  }



  /*//   !!!!!!!!!!!!!!!!!!!!!Cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


  R[0][0] = sqrt(A[0][0]);


  for (i = 1; i < n; i++)
  {
    rr = 0.0;
    rrr = 0.0;
    for (j = 0; j <= i; j++)
    {
      r1 = 0.0;
      r2 = 0.0;
      for (kk = 0; kk <= j - 1; kk++)
      {
        rr = R[i][kk] * R[j][kk] + r1;
        r1 = rr;
      }
      xpivot = R[j][j];
      if (fabs(xpivot) <= 1.e-10)
      {
        printf("nul pivot %d, and R(%d,%d) %.14g\n", j, j, j, xpivot);
        return (*info = 2);
      }
      R[i][j] = (1 / xpivot) * (A[i][j] - rr);
      r3 = 0.0;
      for (kk = 0; kk <= i - 1; kk++)
      {
        rrr = R[i][kk] * R[i][kk] + r3;
        r3 = rrr;
      }
      R[i][i] = sqrt(A[i][i] - rrr);
    }
  }

  /*    //  !!!!!end of cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*   //  !determination of the R tranposeted*/



  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      RT[i][j] = R[j][i];
      invRT[i][j] = 0.;
      invR[i][j] = 0.;
    }
  }


  /*   // !!!!!!!!!inversion of the inf triangular matrix!!!!!!!!!!!!!*/

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (i == j)
      {
        invR[i][j] = 1 / R[i][j];
      }
      else
      {
        invR0 = 0.;
        for (kk = j; kk <= (i - 1); kk++)
        {
          invR[i][j] = (-1 / R[i][i]) * R[i][kk] * invR[kk][j] + invR0;
          invR0 = invR[i][j];
        }
      }
    }
  }

  /*   //  !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*  //  !!!!!!!!!!!!!!!!!!!!!!inversion of the sup triangular matrix!!!!!!!*/

  for (i = 0; i < n; i++)
  {
    invRT[i][i] = 1 / RT[i][i];
  }


  for (i = n - 2; i >= 0; i--)
  {
    for (j = n - 1; j >= 0; j--)
    {
      invRT0 = 0.;
      for (kk = i + 1; kk <= j; kk++)
      {
        invRT[i][j] = (-1 / RT[i][i]) * RT[i][kk] * invRT[kk][j] + invRT0;
        invRT0 = invRT[i][j];
      }
    }
  }

  /*  //  !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/



  /*  //  ! initialisation  s^ = 0*/

  wc = (double*)malloc(n * sizeof(double));
  zc = (double*)malloc(n * sizeof(double));
  znum1 = (double*)malloc(n * sizeof(double));
  wnum1 = (double*)malloc(n * sizeof(double));
  wt = (double*)malloc(n * sizeof(double));
  zt = (double*)malloc(n * sizeof(double));
  kinvnum1 = (double*)malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    wc[i] = 0.0;
    zc[i] = 0.;
    z[i] = 0.;
    w[i] = 0.;
    znum1[i] = 0.;
    wnum1[i] = 0.;
    wt[i] = 0.;
    zt[i] = 0.;
    kinvnum1[i] = 0.;
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      invRTinvR[i][j] = 0.;
      kinv[i][j] = 0.;
    }
  }



  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      invRTinvR0 = 0.;
      for (kk = 0; kk < n; kk++)
      {
        invRTinvR[i][j] = invRT[i][kk] * invR[kk][j] + invRTinvR0;
        invRTinvR0 = invRTinvR[i][j];
      }
    }
  }



  for (i = 0; i < n; i++)
  {
    kinv[i][i] = 1 / k[i][i];
  }


  /*    //  !iteration loops*/

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < *itermax) && (err1 > errmax))
  {

    /*/   !linear stage (zc,wc) -> (z,w)*/

    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, zc, &incx, &beta, wc, &incy);
    dcopy_(&n, qq, &incx, znum1, &incy);
    daxpy_(&n, &alpha, wc, &incx, znum1, &incy);
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, invRTinvR, &n, znum1, &incx, &beta, z, &incy);
    dcopy_(&n, wc, &incx, w, &incy);
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, w, &incy);

    /*    // Local stage (z,w)->(zc,wc)*/


    dcopy_(&n, w, &incx, zt, &incy);
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, zt, &incy);

    for (i = 0; i < n; i++)
    {
      aa = a[i];
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

    dcopy_(&n, wc, &incx, wnum1, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, zt, &incx, wnum1, &incy);
    dcopy_(&n, wnum1, &incx, zt, &incy);

    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, zt, &incx, &beta, zc, &incy);

    /*      // convergence criterium */

    dcopy_(&n, w, &incx, wnum1, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, wc, &incx, wnum1, &incy);
    dcopy_(&n, z, &incx, znum1, &incy);
    daxpy_(&n, &alpha, zc, &incx, znum1, &incy);
    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, znum1, &incx, &beta, wnum1, &incy);

    /*    //    num1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))*/

    num11 = 0.;
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, wnum1, &incx, &beta, kinvnum1, &incy);
    num11 = ddot_(&n, wnum1, &incx, kinvnum1, &incy);

    dcopy_(&n, w, &incx, wnum1, &incy);
    alpha = 1.;
    daxpy_(&n, &alpha, wc, &incx, wnum1, &incy);
    dcopy_(&n, z, &incx, znum1, &incy);
    daxpy_(&n, &alpha, zc, &incx, znum1, &incy);
    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, znum1, &incx, &beta, kinvnum1, &incy);
    den22 = ddot_(&n, znum1, &incx, kinvnum1, &incy);
    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, wnum1, &incx, &beta, kinvnum1, &incy);
    den11 = ddot_(&n, wnum1, &incx, kinvnum1, &incy);

    err0 = num11 / (den11 + den22);
    err1 = sqrt(err0);
    iter1 = iter1 + 1;
    *it_end = iter1;
    *res = err1;

    for (i = 0; i < n; i++)
    {
      /*result_gs[i][iter1-1] = z[i]; */
      fprintf(f101, "%d  %d  %14.7e %14.7e\n", iter1 - 1, i, z[i], w[i]);
    }



  }


  if (err1 > errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter1, err1);
    *info = 1;
  }
  else
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter1, err1);
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
  free(A);
  free(R);
  free(RT);
  free(invRT);
  free(invR);
  free(invRTinvR);
  free(kinv);
  fclose(f101);



}















































