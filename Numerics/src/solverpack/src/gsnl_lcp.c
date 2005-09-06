#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

/*!\file gsnl_lcp.c
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
 * \fn  gsnl_lcp( double *vec , double *q , int *nn , int *itermax , double *tol , double *z , double *omega , int *ispeak ,
 *                double *w , int *it_end , double *res , int *info )
 *
 * gsnl_lcp is a solver for LCP based on the principle of splitting method.
 * (Non Linear Gauss-Seidel)
 *
 * \param vec     Unchanged parameter. Pointer over doubles which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter. Pointer over doubles which contains the components of the right hand side vector.
 * \param nn      Unchanged parameter. Pointer over integer which represents the dimension of the system.
 * \param itermax Unchanged parameter. Pointer over integer which represents the maximum number of iterations allowed.
 * \param tol     Unchanged parameter. Pointer over double which represents the tolerance required.
 * \param omega   Unchanged parameter. Pointer over double which represents the relaxation parameter.
 * \param ispeak  Unchanged parameter. Pointer over integer which represents the output log identifiant
 *                0 - no output
 *                0 < identifiant
 *
 * \param it_end  Modified parameter. Pointer over integer which returns the number of iterations performed by the algorithm.
 * \param res     Modified parameter. Pointer over double which returns the final error value.
 * \param z       Modified parameter. Pointer over doubles which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter. Pointer over doubles which returns the solution of the problem.
 * \param info    Modified parameter. Pointer over integer which returns the termination value
 *                0 - successful
 *                1 - unsuccessful
 *
 * \author Mathieu Renouf
 *
 */

void gsnl_lcp(double *vec , double *q , int *nn , int *itermax , double *tol , double *z , double *omega , int *ispeak ,
              double *w , int *it_end , double *res , int *info)
{


  int n, incx, incy;
  int i, j, iter;

  double qs, err, num, den , zi;
  double *ww, *diag;

  n = *nn;

  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  /* Check for non trivial case */

  incx = 1;
  qs = dnrm2_(&n , &q[0] , &incx);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }
    info = 0;
    return;
  }

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    w[i] = 0.;
  }

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(vec[i * n + i]) < 1e-16)
    {
      info = 2;
      return;
    }
    else diag[i] = 1.0 / vec[i * n + i];
  }

  iter = 0;
  err  = 1.;

  while ((iter < *itermax) && (err > *tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    dcopy_(&n , w , &incx , ww , &incy);

    for (i = 0 ; i < n ; ++i)
    {

      incx = n;
      incy = 1;

      z[i] = 0.0;

      zi = (-q[i] - ddot_(&n , &vec[i] , &incx , z , &incy)) * diag[i];

      if (zi < 0) z[i] = 0.0;
      else z[i] = zi;

      w[i] = (-zi +  z[i]) * vec[i * n + i];

    }

    /* **** Criterium convergence **** */

    qs   = -1;
    incx =  1;
    incy =  1;

    daxpy_(&n , &qs , w , &incx , ww , &incy);

    num = dnrm2_(&n, ww , &incx);
    err = num * den;

    /* **** ********************* **** */

  }

  *it_end = iter;
  *res    = err;


  if (*ispeak > 0)
  {
    if (err > *tol)
    {
      printf(" No convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }

  free(ww);

}
