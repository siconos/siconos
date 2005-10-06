/* Siconos version 1.0, Copyright INRIA 2005.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
/*!\file lcp_nlgs.c
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
 * \fn  lcp_nlgs( int *nn , double *vec , double *q , double *z , double *w , int *info\n,
 *                int *iparamLCP , double *dparamLCP )
 *
 * lcp_nlgs (Non Linear Gauss-Seidel) is a solver for LCP based on the principle of splitting method\n
 *
 * Generic lcp parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param info    Modified parameter which returns the termination value\n
 *                0 - convergence\n
 *                1 - iter = itermax\n
 *                2 - negative diagonal term\n
 *
 * Specific NLGS parameters:\n
 *
 * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
 * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
 *                       0 - no output\n
 *                       0 < active screen output\n
 * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
 *
 * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
 * \param dparamLCP[1] = omega   Input unchanged parameter which represents the relaxation parameter.
 * \param dparamLCP[2] = res     Output modified parameter which returns the final error value.
 *
 * \author Mathieu Renouf
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

void lcp_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info ,
              int *iparamLCP , double *dparamLCP)
{

  int n, incx, incy;
  int i, iter;
  int itermax, ispeak;

  double qs, err, num, den, zi;
  double tol, omega;
  double *ww, *diag;
  double a1, b1;

  char NOTRANS = 'N';

  n = *nn;
  incx = 1;
  incy = 1;

  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  tol   = dparamLCP[0];
  omega = dparamLCP[1];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;

  /* Allocation */

  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  /* Check for non trivial case */

  qs = dnrm2_(&n , q , &incx);

  if (ispeak > 0) printf("\n ||q||= %g \n", qs);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }

    free(ww);
    free(diag);

    *info = 0;
    return;
  }

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    w[i] = 0.;
  }

  /* Intialization of w */

  incx = 1;
  incy = 1;
  dcopy_(&n , q , &incx , w , &incy);

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(vec[i * n + i]) < 1e-16)
    {

      if (ispeak > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem can be solved \n");
      }

      *info = 2;
      free(diag);
      free(ww);

      return;
    }
    else diag[i] = 1.0 / vec[i * n + i];
  }

  /*start iterations*/

  iter = 0;
  err  = 1.;

  incx = 1;
  incy = 1;

  dcopy_(&n , q , &incx , w , &incy);

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    dcopy_(&n , w , &incx , ww , &incy);
    dcopy_(&n , q , &incx , w , &incy);

    for (i = 0 ; i < n ; ++i)
    {

      incx = n;
      incy = 1;

      z[i] = 0.0;

      zi = -(q[i] + ddot_(&n , &vec[i] , &incx , z , &incy)) * diag[i];

      if (zi < 0) z[i] = 0.0;
      else z[i] = zi;

    }

    /* **** Criterium convergence **** */

    incx =  1;
    incy =  1;

    a1 = 1.0;
    b1 = 1.0;

    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , w , &incy);

    qs   = -1.0;
    daxpy_(&n , &qs , w , &incx , ww , &incy);

    num = dnrm2_(&n, ww , &incx);
    err = num * den;

    if (ispeak == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", z[i]);
      for (i = 0 ; i < n ; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* **** ********************* **** */

  }

  iparamLCP[2] = iter;
  dparamLCP[2] = err;

  if (ispeak > 0)
  {
    if (err > tol)
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
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }

  free(ww);
  free(diag);

  return;
}
