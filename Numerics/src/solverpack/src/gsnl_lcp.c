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


/*!\fn int gsnl_lcp(double vec[],double *qq,int *nn,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)


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
\return int : \todo tell whiwh result is good
\author Nineb Sheherazade.
*/
int gsnl_lcp(double vec[], double *qq, int *nn, int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{
  FILE *f13;
  int i, j, kk, iter1, k;
  int n = *nn, incx, incy;
  double errmax = *tol, alpha, beta, mina, maxa, bb, aa;
  double err1, num11, num, den, avn, avt, apn, apt, xn;
  double *q, *zz, *ww, *Mz, *zi, *wi, *y, *yy, *zt, *wnum1;
  char trans;
  double M[n][n];


  q = (double*)malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
    q[i] = qq[i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j];

  yy = (double*)malloc(n * sizeof(double));
  wi = (double*)malloc(n * sizeof(double));
  zi = (double*)malloc(n * sizeof(double));
  ww = (double*)malloc(n * sizeof(double));
  wnum1 = (double*)malloc(n * sizeof(double));
  Mz = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));
  zt = (double*)malloc(n * sizeof(double));


  for (i = 0; i < n; i++)
  {
    y[i] = 0.;
    zi[i] = 0.;
    wi[i] = 0.;
    w[i] = 0.;
    z[i] = 0.;
    Mz[i] = 0.;
    zt[i] = 0.;
    yy[i] = 0.;
    ww[i] = 0.;
  }

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < *itermax) && (err1 > errmax))
  {
    //   !linear stage (zc,wc) -> (z,w)
    iter1 = iter1 + 1;
    incx = 1;
    incy = 1;
    dcopy_(&n, q, &incx, y, &incy);

    trans = 'T';
    alpha = 1.;
    beta = -1.;
    incx = 1;
    incy = 1;
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

      xn = q[i] - avn - apn;

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

    incx = 1;
    incy = 1;
    dcopy_(&n, ww, &incx, y, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, w, &incx, y, &incy);
    dcopy_(&n, y, &incx, wnum1, &incy);
    num = ddot_(&n, wnum1, &incx, wnum1, &incy);
    den = ddot_(&n, ww, &incx, ww, &incy);


    /*  // verif equi
    incx=1;
    incy=1;
    dcopy_(&n,q,&incx,y,&incy);
    alpha = 1.;
    daxpy_(&n,&alpha,w,&incx,y,&incy);
    alpha = 1.;
    beta = -1.;
    dgemv_(&trans,&n,&n,&alpha,M,&n,z,&incx,&beta,y,&incy);

    /*    for (i=0;i<n;i++){
    printf("voici zero (%d) = % g \n",i,y[i]);
    }//*/

    //fin verif equi */


    err1 = sqrt(num) / sqrt(den);
    it_end = &iter1;
    res = &err1;

    printf("iteration numbers %d and error evaluation %e \n ", iter1, err1);
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

  free(Mz);
  free(yy);
  free(wi);
  free(zi);
  free(ww);
  free(wnum1);
  free(y);
  free(q);
  free(zt);
  return *info;
}
