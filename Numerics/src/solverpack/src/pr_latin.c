#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*!\file rp_latin.c


  This subroutine allows the primal resolution of relay problems.

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
-w \in \partial\psi_{[-b, a]}(z)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

*/


double ddot_(int *, double [], int *, double [], int*);


/*!\fn int rp_latin(double vec[],double *qq,int *nn, double * k_latin,double a[],int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   rp_latin is a specific latin solver for primal relay problems.



   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param k_latin On enter a pointer over doubles, the latin coefficient (positive).
   \param a On enter a pointer over doubles, the upper bound.
   \param b On enter a pointer over doubles, the down bound.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param tol On enter a pointer over doubles, the tolerance required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).

  \return int : \todo tell whiwh result is good
   \author Nineb Sheherazade.
 */

rp_latin(double vec[], double *qq, int *nn, double * k_latin, double a[], double b[], int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{

  int i, j, kk, iter1;
  int n = *nn, incx = 1, incy = 1, itt = *itermax;
  double errmax = *tol, alpha, beta, mina;
  double rr, rrr, r1, r2, r3, invR0, invRT0, err1, z0, num11, err0, invRTinvR0;
  double den11, den22;
  double *wc, *zc, *wnum1, *znum1;
  double *zt;
  char trans = 'T';
  //double k[n][n], A[n][n], R[n][n], RT[n][n], invRT[n][n], invR[n][n];
  //double invRTinvR[n][n], kinv[n][n];

  double **k, **A, **R, **RT, **invRT, **invR;
  double  **invRTinvR, **kinv;


  k = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) k[i] = (double*)malloc(n * sizeof(double));
  A = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) A[i] = (double*)malloc(n * sizeof(double));
  R = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) R[i] = (double*)malloc(n * sizeof(double));
  RT = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) RT[i] = (double*)malloc(n * sizeof(double));
  invRT = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) invRT[i] = (double*)malloc(n * sizeof(double));
  invR = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) invR[i] = (double*)malloc(n * sizeof(double));
  invRTinvR = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) invRTinvR[i] = (double*)malloc(n * sizeof(double));
  kinv = (double **)malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) kinv[i] = (double*)malloc(n * sizeof(double));




  printf("itermax %d k_latin %g nn %d vec %g\n", itt, *k_latin, n, vec[0]);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      k[i][j] = 0.;

  for (i = 0; i < n; i++)
    k[i][i] =  *k_latin * vec[i * n + i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      A[i][j] =  vec[i * n + j] + k[i][j];
      R[i][j] = 0.;
    }


  /*  // !!!!!!!!!!!!!!!!!!!!!Cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

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
      if (fabs(R[j][j]) <= 1.e-10)
      {
        printf("nul pivot %d ,and R(%d,%d) %g/n", j, j, j, R[j][j]);
        break;
      }
      R[i][j] = (1 / R[j][j]) * (A[i][j] - rr);
      r3 = 0.0;
      for (kk = 0; kk <= i - 1; kk++)
      {
        rrr = R[i][kk] * R[i][kk] + r3;
        r3 = rrr;
      }
      R[i][i] = sqrt(A[i][i] - rrr);
    }
  }

  /*     // !!!!!end of cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*//  !determination of the R tranposeted*/

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      RT[i][j] = R[j][i];
      invRT[i][j] = 0.;
      invR[i][j] = 0.;
    }


  /*   // !!!!!!!!!inversion of the inf triangular matrix!!!!!!!!!!!!!*/

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (i == j) invR[i][j] = 1 / R[i][j];
      else
      {
        invR0 = 0.;
        for (kk = j; kk <= (i - 1); kk++)
        {
          invR[i][j] = (-1 / R[i][i]) * R[i][kk] * invR[kk][j] + invR0;
          invR0 = invR[i][j];
        }
      }


  /*   // !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*  // !!!!!!!!!!!!!!!!!!!!!!inversion of the sup triangular matrix!!!!!!!*/

  for (i = 0; i < n; i++)
    invRT[i][i] = 1 / RT[i][i];

  for (i = n - 2; i >= 0; i--)
    for (j = n - 1; j >= 0; j--)
    {
      invRT0 = 0.;
      for (kk = i + 1; kk <= j; kk++)
      {
        invRT[i][j] = (-1 / RT[i][i]) * RT[i][kk] * invRT[kk][j] + invRT0;
        invRT0 = invRT[i][j];
      }
    }


  /*// !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


  /*/  ! initialisation  s^ = 0*/

  wc = (double*) malloc(n * sizeof(double));
  zc = (double*) malloc(n * sizeof(double));
  znum1 = (double*) malloc(n * sizeof(double));
  wnum1 = (double*) malloc(n * sizeof(double));
  zt = (double*) malloc(n * sizeof(double));


  for (i = 0; i < n; i++)
  {
    wc[i] = 0.0;
    zc[i] = 0.;
    z[i] = 0.;
    w[i] = 0.;
    znum1[i] = 0.;
    wnum1[i] = 0.;
    zt[i] = 0.;
  }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      invRTinvR[i][j] = 0.;
      kinv[i][j] = 0.;
    }


  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      invRTinvR0 = 0.;
      for (kk = 0; kk < n; kk++)
      {
        invRTinvR[i][j] = invRT[i][kk] * invR[kk][j] + invRTinvR0;
        invRTinvR0 = invRTinvR[i][j];
      }
    }

  for (i = 0; i < n; i++)
    kinv[i][i] = 1 / k[i][i];


  /*    //  !iteration loops*/

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {
    /*    //   !linear stage (zc,wc) -> (z,w)*/

    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, zc, &incx, &beta, wc, &incy);
    dcopy_(&n, qq, &incx, znum1, &incy);
    alpha = 1.;
    daxpy_(&n, &alpha, wc, &incx, znum1, &incy);
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, invRTinvR, &n, znum1, &incx, &beta, z, &incy);
    dcopy_(&n, wc, &incx, w, &incy);
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, w, &incy);

    /*   // Local stage (z,w)->(zc,wc)*/

    dcopy_(&n, z, &incx, zt, &incy);
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, w, &incx, &beta, zt, &incy);

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

    dcopy_(&n, w, &incx, wc, &incy);
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, wc, &incy);
    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, zc, &incx, &beta, wc, &incy);

    /*      // convergence criterium */

    dcopy_(&n, w, &incx, wnum1, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, wc, &incx, wnum1, &incy);
    dcopy_(&n, z, &incx, znum1, &incy);
    daxpy_(&n, &alpha, zc, &incx, znum1, &incy);
    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, znum1, &incx, &beta, wnum1, &incy);
    /*num1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))*/

    num11 = 0.;
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, wnum1, &incx, &beta, znum1, &incy);
    num11 = ddot_(&n, wnum1, &incx, znum1, &incy);

    dcopy_(&n, z, &incx, znum1, &incy);
    daxpy_(&n, &alpha, zc, &incx, znum1, &incy);
    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, znum1, &incx, &beta, wnum1, &incy);
    den22 = ddot_(&n, znum1, &incx, wnum1, &incy);
    dcopy_(&n, w, &incx, wnum1, &incy);
    alpha = 1.;
    daxpy_(&n, &alpha, wc, &incx, wnum1, &incy);
    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, wnum1, &incx, &beta, znum1, &incy);
    den11 = ddot_(&n, wnum1, &incx, znum1, &incy);
    err0 = num11 / (den11 + den22);
    err1 = sqrt(err0);
    iter1 = iter1 + 1;
    *it_end = iter1;
    *res = err1;
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
  free(zt);

  for (i = 0; i < n; i++) free(k[i]);
  free(k);
  for (i = 0; i < n; i++) free(A[i]);
  free(A);
  for (i = 0; i < n; i++) free(R[i]);
  free(R);
  for (i = 0; i < n; i++) free(RT[i]);
  free(RT);
  for (i = 0; i < n; i++) free(invRT[i]);
  free(invRT);
  for (i = 0; i < n; i++) free(invR[i]);
  free(invR);
  for (i = 0; i < n; i++) free(invRTinvR[i]);
  free(invRTinvR);
  for (i = 0; i < n; i++) free(kinv[i]);
  free(k);



}
