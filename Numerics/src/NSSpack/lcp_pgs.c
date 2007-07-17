/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2006.
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
/*!\file lcp_pgs.c

  This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
  Try \f$(z,w)\f$ such that:\n
  \f$
   \left\lbrace
   \begin{array}{l}
     w - M z = q\\
     0 \le z \perp w \ge 0\\
    \end{array}
   \right.
  \f$

  where M is an (\f$nn \times nn\f$)-matrix, q , w and z nn-vectors.
*/
/*!\fn  void lcp_pgs( int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP )

  lcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for LCP.\n

  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
                 0 : convergence\n
                 1 : iter = itermax\n
                 2 : negative diagonal term

  \param iparamLCP  On enter/return a vector of integers:\n
                - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
                - iparamLCP[1] = verbose  On enter, the output log identifiant:\n
                        0 : no output\n
                        >0: active screen output\n
                - iparamLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.

  \param dparamLCP  On enter/return a vector of doubles:\n
                - dparamLCP[0] = tol     On enter, the tolerance required.
                - dparamLCP[1] = omega   On enter, the relaxation parameter (not yet available).
                - dparamLCP[2] = res     On return, the final error value.

  \author Mathieu Renouf

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blaslapack.h"
#include <math.h>
int lcp_compute_error(int n, double *vec , double *q , double *z , int chat, double *w, double *err);

void lcp_pgs(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{


  int n, incx, incy;
  int i, iter;
  int itermax, verbose;
  int incxn;
  double qs, err, num, den, zi;
  double tol, omega;
  double *ww, *diag;
  double a1, b1;

  char NOTRANS = 'N';

  n = *nn;
  incx = 1;
  incy = 1;
  incxn = n;
  /* Recup input */

  itermax = iparamLCP[0];
  verbose  = iparamLCP[1];

  tol   = dparamLCP[0];
  omega = dparamLCP[1];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;

  /* Allocation */

  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  /* Check for non trivial case */
  integer totoincx;
  totoincx = (integer)incx;

  qs = dnrm2_((integer *)&n , q , (integer *)&incx);

  if (verbose > 0) printf("\n ||q||= %g \n", qs);

  den = 1.0 / qs;

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    w[i] = 0.;
  }

  /* Intialization of w */

  incx = 1;
  incy = 1;
  dcopy_((integer *)&n , q , (integer *)&incx , w , (integer *)&incy);

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(vec[i * n + i]) < 1e-16)
    {

      if (verbose > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
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

  dcopy_((integer *)&n , q , (integer *)&incx , w , (integer *)&incy);

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    dcopy_((integer *)&n , w , (integer *)&incx , ww , (integer *)&incy);
    dcopy_((integer *)&n , q , (integer *)&incx , w , (integer *)&incy);

    for (i = 0 ; i < n ; ++i)
    {

      incx = n;
      incy = 1;

      z[i] = 0.0;

      zi = -(q[i] + ddot_((integer *)&n , &vec[i] , (integer *)&incx , z , (integer *)&incy)) * diag[i];

      if (zi < 0) z[i] = 0.0;
      else z[i] = zi;

      /* z[i]=fmax(0.0,-( q[i] + ddot_( (integer *)&n , &vec[i] , (integer *)&incxn , z , (integer *)&incy ))*diag[i]);*/

    }

    /* **** Criterium convergence **** */

    /*     incx =  1; */
    /*     incy =  1; */

    /*     a1 = 1.0; */
    /*     b1 = 1.0; */

    /*     dgemv_( &NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , (integer *)&incx , &b1 , w , (integer *)&incy ); */

    /*     qs   = -1.0; */
    /*     daxpy_( (integer *)&n , &qs , w , (integer *)&incx , ww , (integer *)&incy ); */

    /*     num = dnrm2_( (integer *)&n, ww , (integer *)&incx ); */

    /*     err = num*den; */

    /* **** Criterium convergence compliant with filter_result_LCP **** */



    lcp_compute_error(n, vec, q, z, verbose, w, &err);

    //err = err ;



    if (verbose == 2)
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

  if (err > tol)
  {
    printf("Siconos/Numerics: lcp_pgs: No convergence of PGS after %d iterations\n" , iter);
    printf("Siconos/Numerics: lcp_pgs: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: lcp_pgs: Convergence of PGS after %d iterations\n" , iter);
      printf("Siconos/Numerics: lcp_pgs: The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(ww);
  free(diag);

  return;
}
