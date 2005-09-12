#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

/*!\file lcp_cpg.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).
 * Try \f$(z,w)\f$ such that:
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    M z + q= w\\
 *    0 \le z \perp w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * where M is an (n x n)-matrix, q , w and z n-vectors.
 *
 *!\fn  lcp_cpg( double *vec , double *q , int *nn , int *itermax , double *tol , int *ispeak , double *z ,
 *               double *w , int *it_end , double *res , int *info)
 *
 * lcp_cpg is a basic cpg (conjugated projected gradient) solver for LCP.
 *
 * \param double* vec    Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param double* q      Unchanged parameter which contains the components of the right hand side vector.
 * \param int *nn        Unchanged parameter which represents the dimension of the system.
 * \param int *itermax   Unchanged parameter which represents the maximum number of iterations allowed.
 * \param double *tol    Unchanged parameter which represents the tolerance required.
 * \param double *omega  Unchanged parameter which represents the relaxation parameter.
 * \param int *ispeak    Unchanged parameter which represents the output log identifiant
 *                       0 - no output
 *                       0 < identifiant
 *
 * \param int* it_end    Modified parameter which returns the number of iterations performed by the algorithm.
 * \param double* res    Modified parameter which returns the final error value.
 * \param double* z      Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param double* w      Modified parameter which returns the solution of the problem.
 * \param int* info      Modified parameter which returns the termination value
 *                       0 - convergence
 *                       1 - iter = itermax
 *                       2 - negative diagonal term
 *                       3 - pWp nul
 *
 * \author Mathieu Renouf
 *
 */

void lcp_cpg(double *vec , double *q , int *nn , int *itermax , double *tol , int *ispeak , double *z ,
             double *w , int *it_end , double * res , int *info)
{


  int n, incx, incy;
  int i, iter;

  double err, a1, b1 , qs;

  double alpha, beta, rp, pMp;
  double den, num;

  char NOTRANS = 'N';

  int *status;
  double *zz , *dz , *pp , *rr, *ww, *Mp;

  status = (int*)malloc(n * sizeof(int));

  dz = (double*)malloc(n * sizeof(double));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  *info = 1;
  n = *nn;
  incx = 1;

  qs = dnrm2_(&n , &q[0] , &incx);

  //printf( " Norm: %g \n", qs );

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }
    *info = 0;
    return;
  }

  for (i = 0 ; i < n ; ++i) ww[i] = 0.;

  for (i = 0; i < n; ++i)
  {

    status[i] = 0;

    ww[i] = 0.;
    rr[i] = 0.;
    pp[i] = 0.;
    zz[i] = 0.;

    Mp[i] = 0.;

  }

  /* rr = -Wz + q */
  incx = 1;
  incy = 1;

  dcopy_(&n , q , &incx , rr , &incy);

  a1 = -1.;
  b1 = -1.;

  dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , rr , &incy);

  /* Initialization of gradients */
  /* rr -> p and rr -> w */

  dcopy_(&n , rr , &incx , ww , &incy);
  dcopy_(&n , rr , &incx , pp , &incy);

  iter = 0.0;
  err  = 1.0 ;

  while ((iter < *itermax) && (err > *tol))
  {

    ++iter;

    /* Compute initial pMp */

    incx = 1;
    incy = 1;

    dcopy_(&n , pp , &incx , Mp , &incy);

    a1 = 1.0;
    b1 = 0.0;

    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , Mp , &incx , &b1 , w , &incy);

    pMp = ddot_(&n , pp , &incx , w  , &incy);

    if (fabs(pMp) < 1e-16)
    {

      if (*ispeak > 0)
      {
        printf(" Operation no conform at the iteration %d \n", iter);
        printf(" Alpha can be obtained with pWp = %10.4g  \n", pMp);
      }

      *it_end = iter;
      *res    = err;

      return (*info = 3);
    }

    rp  = ddot_(&n , pp , &incx , rr , &incy);

    alpha = rp / pMp;

    /*
     * Iterate prediction:
     * z' = z + alpha*p
     *
     */

    dcopy_(&n , z , &incx , dz , &incy);

    daxpy_(&n , &alpha , pp , &incx , z , &incy);

    /* Iterate projection*/

    for (i = 0; i < n; ++i)
    {
      if (z[i] > 0.0)
      {
        status[i] = 1;
      }
      else
      {
        z[i] = 0.0;
        status[i] = 0;
      }
    }

    /* rr = -Wz + q */

    dcopy_(&n , q , &incx , rr , &incy);

    a1 = -1.;
    b1 = -1.;

    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , rr , &incy);

    /* Gradients projection
     * rr --> ww
     * pp --> zz
     */

    for (i = 0; i < n; ++i)
    {

      if (status[i])
      {
        ww[i] = rr[i];
        zz[i] = pp[i];
      }
      else
      {
        if (rr[i] < 0)
        {
          ww[i] = 0.0;
          zz[i] = 0.0;
        }
        else
        {
          ww[i] = rr[i];
          if (pp[i] < 0) zz[i] = 0.0;
          else zz[i] = pp[i];
        }
      }
    }

    /*   beta = -w.Mp / pMp  */

    rp = ddot_(&n , ww , &incx, w , &incy);

    beta = -rp / pMp;

    dcopy_(&n , ww , &incx , pp , &incy);
    daxpy_(&n, &beta , zz , &incx , pp , &incy);

    a1 = -1;
    daxpy_(&n , &a1 , z , &incx , dz , &incy);
    num = dnrm2_(&n , dz , &incx);
    err = num * den;

  }

  *it_end = iter;
  *res    = err;

  if (*ispeak > 0)
  {
    if (err > *tol)
    {
      printf(" No convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else if (err < *tol) *info = 0;

  free(Mp);

  free(ww);
  free(rr);
  free(pp);
  free(zz);

  free(dz);

}
