#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*!\file gsnl_lcp.c


This subroutine allows the resolution of LCP (Linear Complementary Problem).
Try \f$(z,w)\f$ such that:

\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z \perp w \ge 0\\
\end{array}
\right.
\f$

here M is an n by n  matrix, q an n-dimensional vector, w an n-dimensional  vector and z an n-dimensional vector.
*/


double ddot_(int *, double [], int *, double [], int*);


/*!\fn  gsnl_lcp(double vec[],double *qq,int *nn,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)


gsnl_lcp is a basic gsnl (Gauss-Seidel Non Linear) solver for LCP.

\param vec On enter a pointer over doubles containing the components of the double matrix with a fortran90 allocation.
\param qq On enter a pointer over doubles containing the components of the double vector.
\param nn On enter a pointer over integers, the dimension of the second member.
\param itermax On enter a pointer over integers, the maximum iterations required.
\param tol On enter a pointer over doubles, the tolerance required.
\param it_end On enter a pointer over integers, the number of iterations carried out.
\param res On return a pointer over doubles, the error value.
\param z On return double vector, the solution of the problem.
\param w On return double vector, the solution of the problem.
\param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).

\author Nineb Sheherazade.
*/



gsnl_lcp(double vec[], double *qq, int *nn, int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{
  int i, j, kk, iter1, k, itt = *itermax;
  int n = *nn, incx = 1, incy = 1;
  double errmax = *tol, alpha, beta;
  double err1, num11, num, den, avn, avt, apn, apt, xn;
  double *q, *zz, *ww, *y;
  char trans = 'T';
  double M[n][n];



  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j]; /* vec has a F90 storage and M has a C storage*/



  ww = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));



  for (i = 0; i < n; i++)
  {
    y[i] = 0.;
    w[i] = 0.;
    z[i] = 0.;
    ww[i] = 0.;
  }

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {
    iter1 = iter1 + 1;
    dcopy_(&n, qq, &incx, y, &incy);


    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, ww, &incy);

    for (i = 0; i < n; i++)
    {
      avn = 0.;
      apn = 0.;

      for (j = 0; j <= i - 1; j++)
        avn = avn + M[i][j] * z[j];

      for (k = i + 1; k < n; k++)
        apn = apn + M[i][k] * z[k];

      xn = qq[i] - avn - apn;

      if (xn > 0.0)
      {
        z[i] = xn / M[i][i];
        w[i] = 0.;
      }
      else
      {
        z[i] = 0.;
        w[i] = -xn;
      }
    }


    /* ///////// Criterium convergence ///////////// */


    alpha = -1.;
    daxpy_(&n, &alpha, w, &incx, ww, &incy);
    num = ddot_(&n, ww, &incx, ww, &incy);
    den = ddot_(&n, qq, &incx, qq, &incy);

    err1 = sqrt(num) / sqrt(den);
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



  free(ww);
  free(y);
}
