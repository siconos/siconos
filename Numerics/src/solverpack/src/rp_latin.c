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
-w \in \partial\psi_{[-a, a]}(z)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

*/


//extern ddotf_(int *,double [],double [],double*);

double ddot_(int *, double [], int *, double [], int*);


/*!\fn int rp_latin(double vec[],double *qq,int *nn, double * k_latin,double a[],int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   rp_latin is a specific latin solver for primal relay problems.



   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param k_latin On enter a pointer over doubles, the latin coefficient (positive).
   \param a On enter a pointer over doubles, the bound.
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
int rp_latin(double vec[], double *qq, int *nn, double * k_latin, double a[], int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{
  FILE *f13;
  int i, j, kk, iter1;
  int n = *nn, incx, incy, info3;
  double errmax = *tol, alpha, beta, mina, maxa, bb, cc, zw, aa;
  double rr, rrr, r1, r2, r3, invR0, invRT0, err1, z0, num11, err0, invRTinvR0;
  double den11, den22, vv;
  double *q, *resveclat, *wc, *zc, *kinvden1, *kinvden2, *wt, *maxwt, *y, *wnum1, *znum1;
  double *zt, *maxzt, *num1, *kinvnum1, *den1, *den2, *wden1, *zden1, *Rchole;
  char trans;
  double k[n][n], Mtp[n][n], M[n][n], A[n][n], R[n][n], RT[n][n], invRT[n][n], invR[n][n];
  double invRTinvR[n][n], kinv[n][n], xx[n][*itermax], kinvwden1[n], kzden1[n];


  q = (double*)malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
    q[i] = -qq[i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      M[i][j] = vec[i * n + j];
      Mtp[i][j] = 0.;
    }

  for (i = 0; i < n; i++)
    Mtp[i][i] = M[i][i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      k[i][j] = 0.;

  for (i = 0; i < n; i++)
    k[i][i] = *k_latin * Mtp[i][i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = M[i][j] + k[i][j];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      R[i][j] = 0.;

  // !!!!!!!!!!!!!!!!!!!!!Cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  // !!!!!end of cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //  !determination of the R tranposeted
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      RT[i][j] = R[j][i];

  // !inverse of R and RT
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      invRT[i][j] = 0.;
      invR[i][j] = 0.;
    }

  // !!!!!!!!!inversion of the inf triangular matrix!!!!!!!!!!!!!
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
  // !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // !!!!!!!!!!!!!!!!!!!!!!inversion of the sup triangular matrix!!!!!!!
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
  // !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  resveclat = (double*)malloc(*itermax * sizeof(double));

  for (i = 0; i < n; i++)
    resveclat[i] = 0.;

  //  ! initialisation  s^ = 0
  wc = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));
  zc = (double*)malloc(n * sizeof(double));
  znum1 = (double*)malloc(n * sizeof(double));
  wnum1 = (double*)malloc(n * sizeof(double));
  kinvden1 = (double*)malloc(n * sizeof(double));
  kinvden2 = (double*)malloc(n * sizeof(double));
  wt = (double*)malloc(n * sizeof(double));
  maxwt = (double*)malloc(n * sizeof(double));
  zt = (double*)malloc(n * sizeof(double));
  maxzt = (double*)malloc(n * sizeof(double));
  num1 = (double*)malloc(n * sizeof(double));
  kinvnum1 = (double*)malloc(n * sizeof(double));
  den1 = (double*)malloc(n * sizeof(double));
  den2 = (double*)malloc(n * sizeof(double));
  wden1 = (double*)malloc(n * sizeof(double));
  zden1 = (double*)malloc(n * sizeof(double));


  for (i = 0; i < n; i++)
  {
    resveclat[i] = 0.;
    wc[i] = 0.0;
    y[i] = 0.;
    zc[i] = 0.;
    z[i] = 0.;
    w[i] = 0.;
    znum1[i] = 0.;
    wnum1[i] = 0.;
    kinvden1[i] = 0.;
    kinvden2[i] = 0.;
    wt[i] = 0.;
    maxzt[i] = 0.;
    maxwt[i] = 0.;
    zt[i] = 0.;
    num1[i] = 0.;
    kinvnum1[i] = 0.;
    den1[i] = 0.;
    den2[i] = 0.;
  }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      invRTinvR[i][j] = 0.;
      kinv[i][j] = 0.;
    }

  //  !iteration loops
  iter1 = 0;
  err1 = 1.;

  while ((iter1 < *itermax) && (err1 > errmax))
  {
    //   !linear stage (zc,wc) -> (z,w)
    incx = 1;
    incy = 1;
    dcopy_(&n, wc, &incx, y, &incy);

    trans = 'T';
    alpha = 1.;
    beta = 1.;
    incx = 1;
    incy = 1;

    dgemv_(&trans, &n, &n, &alpha, k, &n, zc, &incx, &beta, y, &incy);

    incx = 1;
    incy = 1;
    dcopy_(&n, y, &incx, wt, &incy);
    //    wt(:) = wc(:) + matmul(k,zc)

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
    {
      z0 = 0.;
      vv = 0.;
      for (j = 0; j < n; j++)
      {
        vv = wt[j] - q[j];
        z[i] = invRTinvR[i][j] * vv + z0;
        z0 = z[i];
      }
    }


    incx = 1;
    incy = 1;
    dcopy_(&n, qq, &incx, y, &incy);
    alpha = 1.;
    daxpy_(&n, &alpha, wt, &incx, y, &incy);
    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, A, &n, z, &incx, &beta, y, &incy);



    for (i = 0; i < n; i++)
      xx[i][iter1] = z[i];

    incx = 1;
    incy = 1;
    dcopy_(&n, wt, &incx, y, &incy);

    trans = 'T';
    alpha = -1.;
    beta = 1.;
    incx = 1;
    incy = 1;

    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, w, &incy);

    //   w(:) = wt(:)-matmul(k(:,:),z(:))



    // Local stage (z,w)->(zc,wc)

    for (i = 0; i < n; i++)
      kinv[i][i] = 1 / k[i][i];

    dcopy_(&n, z, &incx, y, &incy);
    alpha = -1.;
    beta = 1.;
    incx = 1;
    incy = 1;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, w, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, zt, &incy);

    for (i = 0; i < n; i++)
    {
      aa = a[i];
      minf(&aa, &zt[i], &mina);
      bb = -aa;
      maxf(&mina, &bb, &maxa);
      zc[i] = maxa;
    }


    dcopy_(&n, w, &incx, y, &incy);

    trans = 'T';
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, y, &incy);


    incx = 1;
    incy = 1;
    dcopy_(&n, y, &incx, wt, &incy);
    dcopy_(&n, wt, &incx, y, &incy);

    trans = 'T';
    alpha = 1.;
    beta = 1.;
    incx = 1;
    incy = 1;


    dgemv_(&trans, &n, &n, &alpha, k, &n, zc, &incx, &beta, y, &incy);


    incx = 1;
    incy = 1;
    dcopy_(&n, y, &incx, wc, &incy);


    // convergence criterium
    incx = 1;
    incy = 1;
    dcopy_(&n, w, &incx, y, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, wc, &incx, y, &incy);
    dcopy_(&n, y, &incx, wnum1, &incy);
    dcopy_(&n, z, &incx, y, &incy);
    daxpy_(&n, &alpha, zc, &incx, y, &incy);
    dcopy_(&n, y, &incx, znum1, &incy);

    dcopy_(&n, wnum1, &incx, y, &incy);
    trans = 'T';
    alpha = 1.;
    beta = 1.;


    dgemv_(&trans, &n, &n, &alpha, k, &n, znum1, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, num1, &incy);

    //    num1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))

    num11 = 0.;
    incx = 1;
    incy = 1;
    dcopy_(&n, num1, &incx, y, &incy);
    trans = 'T';
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, num1, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, kinvnum1, &incy);

    num11 = ddot_(&n, num1, &incx, kinvnum1, &incy);

    // rectif ici
    incx = 1;
    incy = 1;
    dcopy_(&n, w, &incx, y, &incy);
    // trans='t';
    alpha = 1.;
    daxpy_(&n, &alpha, wc, &incx, y, &incy);
    dcopy_(&n, y, &incx, wden1, &incy);
    dcopy_(&n, z, &incx, y, &incy);
    // trans='t';
    daxpy_(&n, &alpha, zc, &incx, y, &incy);

    dcopy_(&n, y, &incx, zden1, &incy);

    dcopy_(&n, wden1, &incx, y, &incy);
    //   trans='t';
    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, zden1, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, kzden1, &incy);

    den22 = ddot_(&n, zden1, &incx, kzden1, &incy);

    dcopy_(&n, wden1, &incx, y, &incy);
    //   trans='t';
    beta = 0.;
    alpha = 1.;


    dgemv_(&trans, &n, &n, &alpha, kinv, &n, wden1, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, kinvwden1, &incy);

    den11 = ddot_(&n, wden1, &incx, kinvwden1, &incy);

    //// rectif ici
    err0 = num11 / (den11 + den22);
    err1 = sqrt(err0);
    resveclat[iter1] = err1;
    iter1 = iter1 + 1;
    it_end = &iter1;
    res = &err1;
  }

  if (err1 > errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter1, err1);
  }

  if (err1 < errmax)
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter1, err1);
  }

  if (err1 < errmax) *info = 0;
  else *info = 1;

  free(q);
  free(resveclat);
  free(wc);
  free(y);
  free(zc);
  free(znum1);
  free(wnum1);
  free(kinvden1);
  free(kinvden2);
  free(wt);
  free(maxwt);
  free(zt);
  free(maxzt);
  free(num1);
  free(kinvnum1);
  free(den1);
  free(den2);
  free(wden1);
  free(zden1);
  return *info ;
}
