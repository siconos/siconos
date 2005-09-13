#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

/*!\file lcp_cpg.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
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
 * \fn  lcp_cpg( int *nn , double *vec , double *q , double *z , int *info ,
 *               int *iparamLCP , double *dparamLCP )
 *
 * lcp_cpg is a cpg (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.\n
 *
 * Generic lcp parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param info    Modified parameter which returns the termination value\n
 *                0 - convergence
 *                1 - iter = itermax
 *                2 - negative diagonal term
 *                3 - pWp nul
 *
 * Specific CPG parameters:\n
 *
 * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
 * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
 *                       0 - no output\n
 *                       0 < active screen ouput\n
 * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
 *
 * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
 * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
 *
 * \author Mathieu Renouf
 *
 */

void lcp_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info ,
             int *iparamLCP , double *dparamLCP)
{

  int n, incx, incy;
  int i, iter;
  int itermax, ispeak;

  double err, a1, b1 , qs;

  double alpha, beta, rp, pMp;
  double den, num, tol;

  char NOTRANS = 'N';

  int *status;
  double *zz , *dz , *pp , *rr, *ww, *Mp;

  *info = 1;
  n = *nn;

  /*input*/

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  tol   = dparamLCP[0];

  /*output*/

  iparamLCP[2] = 0;
  dparamLCP[1] = 0.0;

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

  /* Allocations */

  status = (int*)malloc(n * sizeof(int));

  dz = (double*)malloc(n * sizeof(double));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  incx = 1;

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

  while ((iter < itermax) && (err > tol))
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

      if (ispeak > 0)
      {
        printf(" Operation no conform at the iteration %d \n", iter);
        printf(" Alpha can be obtained with pWp = %10.4g  \n", pMp);
      }

      free(Mp);
      free(ww);
      free(rr);
      free(pp);
      free(zz);
      free(dz);

      iparamLCP[2] = iter;
      dparamLCP[1] = err;

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

  iparamLCP[2] = iter;
  dparamLCP[1] = err;

  if (ispeak > 0)
  {
    if (err > tol)
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
  else if (err < tol) *info = 0;

  free(Mp);

  free(ww);
  free(rr);
  free(pp);
  free(zz);

  free(dz);

}
