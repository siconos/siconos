#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*!\file gcp_lcp.c


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


/*!\fn int gcp_lcp(double vec[],double *qq,int *nn,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   gcp_lcp is a basic gcp (gradient conjugated projected) solver for LCP.

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
int gcp_lcp(double vec[], double *qq, int *nn, int * itermax, double * tol, double z[], double wsol[], int *it_end, double * res, int *info)
{
  FILE *f13;
  int i, j, kk, iter1, k;
  int n = *nn, incx, incy;
  double errmax = *tol, alpha, beta, mina, maxa, bb, aa, rp;
  double err1, num11, num, den, avn, avt, apn, apt, xn, pmp, alpha1, beta1;
  double *q, *zz, *ww, *Mz, *zi, *wi, *y, *p, *zt, *wnum1, *r, *w, *Mp, *v;
  char trans;
  double M[n][n];



  q = (double*)malloc(n * sizeof(double));
  w = (double*)malloc(n * sizeof(double));
  r = (double*)malloc(n * sizeof(double));
  v = (double*)malloc(n * sizeof(double));



  for (i = 0; i < n; i++)
    q[i] = qq[i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j];

  p = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));
  wi = (double*)malloc(n * sizeof(double));
  zi = (double*)malloc(n * sizeof(double));
  ww = (double*)malloc(n * sizeof(double));
  wnum1 = (double*)malloc(n * sizeof(double));
  Mz = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));
  zt = (double*)malloc(n * sizeof(double));
  Mp = (double*)malloc(n * sizeof(double));



  for (i = 0; i < n; i++)
  {
    y[i] = 0.;
    zz[i] = 0.;
    zi[i] = 0.;
    wi[i] = 0.;
    w[i] = 0.;
    z[i] = 0.;
    r[i] = 0.;
    wsol[i] = 0.;
    v[i] = 0.;
    Mz[i] = 0.;
    zt[i] = 0.;
    p[i] = 0.;
    Mp[i] = 0.;
    ww[i] = 0.;
  }

  incx = 1;
  incy = 1;
  dcopy_(&n, q, &incx, y, &incy);

  trans = 'T';
  alpha = -1.;
  beta = 1.;
  dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy);
  dcopy_(&n, y, &incx, r, &incy);


  incx = 1;
  incy = 1;
  dcopy_(&n, r, &incx, w, &incy);


  incx = 1;
  incy = 1;
  dcopy_(&n, r, &incx, p, &incy);

  trans = 'T';
  alpha = 1.;
  beta = 0.;
  dgemv_(&trans, &n, &n, &alpha, M, &n, p, &incx, &beta, y, &incy);
  dcopy_(&n, y, &incx, Mp, &incy);

  pmp = ddot_(&n, p, &incx, Mp, &incy);

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < *itermax) && (err1 > errmax))
  {
    incx = 1;
    incy = 1;
    dcopy_(&n, r, &incx, v, &incy);

    iter1 = iter1 + 1;
    if (pmp == 0.)
    {
      printf("operation other alpha not conform at the iteration %d", iter1);
      return (*info = 3) ;
    }

    // alpha1 = r.p / pmp
    rp = ddot_(&n, p, &incx, r, &incy);
    alpha1 = rp / pmp;

    incx = 1;
    incy = 1;
    dcopy_(&n, z, &incx, y, &incy);

    alpha = alpha1;
    daxpy_(&n, &alpha, p, &incx, y, &incy);
    dcopy_(&n, y, &incx, zi, &incy);


    for (j = 0; j < n; j++)
      if (zi[j] <= 0.) z[j] = 0.;
      else z[j] = zi[j];

    incx = 1;
    incy = 1;
    dcopy_(&n, q, &incx, y, &incy);

    trans = 'T';
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, r, &incy);


    for (j = 0; j < n; j++)
    {
      if (z[j] <= 0.)
      {
        if (r[j] <= 0.) w[j] = 0.;
        else w[j] = r[j];
      }
      else
      {
        w[j] = r[j];
      }
    }



    for (j = 0; j < n; j++)
    {
      if (z[j] <= 0.)
      {
        if (p[j] <= 0.) zz[j] = 0.;
        else zz[j] = p[j];
      }
      else zz[j] = p[j];
    }

    // beta = -w.Mp / pmp

    rp = ddot_(&n, w, &incx, Mp, &incy);
    beta1 = -rp / pmp;


    dcopy_(&n, w, &incx, y, &incy);
    alpha = beta1;
    daxpy_(&n, &alpha, zz, &incx, y, &incy);
    dcopy_(&n, y, &incx, p, &incy);


    trans = 'T';
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, p, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, Mp, &incy);

    pmp = ddot_(&n, p, &incx, Mp, &incy);


    /* ///////// Criterium convergenc ///////////// */
    // normr=sqrt(dot_product(r(:)-v(:),r(:)-v(:)))/sqrt(dot_product(v(:),v(:)))

    incx = 1;
    incy = 1;
    dcopy_(&n, r, &incx, y, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, v, &incx, y, &incy);
    dcopy_(&n, y, &incx, wnum1, &incy);
    num = ddot_(&n, wnum1, &incx, wnum1, &incy);
    den = ddot_(&n, v, &incx, v, &incy);
    err1 = sqrt(num) / sqrt(den);

    it_end = &iter1;
    res = &err1;

    alpha = -1.;
    daxpy_(&n, &alpha, r, &incx, y, &incy);
    dcopy_(&n, y, &incx, wsol, &incy);

    printf("iteration numbers %d and error evaluation %e \n ", iter1, err1);

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
  }

  free(q);
  free(w);
  free(r);
  free(v);
  free(Mp);
  free(p);
  free(zz);
  free(wi);
  free(zi);
  free(ww);
  free(wnum1);
  free(Mz);
  free(y);
  free(zt);
  return *info;
}
