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


/*!\fn  gcp_lcp(double vec[],double *q,int *nn,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   gcp_lcp is a basic gcp (gradient conjugated projected) solver for LCP.

   \param vec On enter a pointer over doubles containing the components of the double matrix with a fortran90 allocation.
   \param q On enter a pointer over doubles containing the components of the double vector.
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





gcp_lcp(double vec[], double *q, int *nn, int * itermax, double * tol, double z[], double *vv, int *it_end, double * res, int *info)
{
  int i, j, iter1;
  int n = *nn, incx = 1, incy = 1, itt = *itermax;
  double errmax = *tol, alpha, beta, rp;
  double err1, pMp, alpha1, beta1, comp;
  double *zz, *y, *p, *r, *w, *Mp;
  char trans = 'T';
  double M[n][n];






  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j];



  w = (double*)malloc(n * sizeof(double));
  r = (double*)malloc(n * sizeof(double));
  p = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));
  Mp = (double*)malloc(n * sizeof(double));



  for (i = 0; i < n; i++)
  {
    y[i] = 0.;
    zz[i] = 0.;
    w[i] = 0.;
    z[i] = 0.;
    r[i] = 0.;
    vv[i] = 0.;
    p[i] = 0.;
    Mp[i] = 0.;
  }


  dcopy_(&n, q, &incx, y, &incy);


  alpha = -1.;
  beta = 1.;
  dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy); /*//r*/
  dcopy_(&n, y, &incx, r, &incy);



  dcopy_(&n, r, &incx, p, &incy);


  alpha = 1.;
  beta = 0.;
  dgemv_(&trans, &n, &n, &alpha, M, &n, p, &incx, &beta, Mp, &incy);

  pMp = ddot_(&n, p, &incx, Mp, &incy);

  rp = ddot_(&n, p, &incx, r, &incy);

  iter1 = 0;
  err1 = pMp;

  while ((iter1 < itt) && (err1 >= errmax))
  {

    iter1 = iter1 + 1;

    if (pMp == 0.)
    {
      printf("operation other alpha not conform at the iteration %d", iter1);
      return (*info = 3) ;
    }


    alpha1 = rp / pMp;



    alpha = alpha1;
    daxpy_(&n, &alpha, p, &incx, z, &incy);


    for (j = 0; j < n; j++)
      if (z[j] <= 0.) z[j] = 0.;


    dcopy_(&n, q, &incx, y, &incy);

    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, r, &incy);
    dcopy_(&n, r, &incx, w, &incy);
    dcopy_(&n, p, &incx, zz, &incy);

    for (j = 0; j < n; j++)
    {
      if (z[j] == 0. && w[j] < 0)
      {
        w[j] = 0.;
        zz[j] = 0.;
      }
    }


    /*   beta = -w.Mp / pMp  */
    rp = ddot_(&n, w, &incx, Mp, &incy);
    beta1 = -rp / pMp;



    alpha = beta1;
    daxpy_(&n, &alpha, zz, &incx, w, &incy);
    dcopy_(&n, w, &incx, p, &incy);



    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, p, &incx, &beta, Mp, &incy);


    pMp = ddot_(&n, p, &incx, Mp, &incy);


    /* ///////// Criterium convergence : err1= pMp ///////////// */



    err1 = pMp;


    dcopy_(&n, q, &incx, y, &incy);
    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, vv, &incy);
    rp = ddot_(&n, p, &incx, r, &incy);

    *it_end = iter1;
    *res = err1;


    /*    printf("iteration numbers %d and error evaluation %g \n ",iter1,err1);   */
  }


  if (err1 >= errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter1, err1);
    *info = 1;
  }
  else
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter1, err1);
    comp = ddot_(&n, z, &incx, vv, &incy);
    printf(" the complementarity %g\n", comp);
    *info = 0;
  }



  free(w);
  free(r);
  free(Mp);
  free(p);
  free(zz);
  free(y);


}
