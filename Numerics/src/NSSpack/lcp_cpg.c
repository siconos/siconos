/* Siconos-Numerics version 2.1.0, Copyright INRIA 2005-2006.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

/*!\file lcp_cpg.c

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
/*!\fn  void lcp_cpg( int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP )

  lcp_cpg is a CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.


  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which contains the termination value:
                  0: convergence\n
                  1: iter = itermax\n
                  2: negative diagonal term\n
                  3: pWp nul

 \param iparamLCP On enter/return, a vector of integers:
                - iparamLCP[0] = itermax On enter, an integer which represents the maximum number of iterations allowed.
                - iparamLCP[1] = ispeak  On enter, an integer which represents the output log identifiant:
                   0 : no output\n
                   > 0: active screen ouput
                - iparamLCP[2] = it_end  On return, a double which represents the number of iterations performed by the algorithm.

 \peram dparamLCP On enter/return, a vector of doubles:
                - dparamLCP[0] = tol     On enter, a double which represents the tolerance required.
                - dparamLCP[1] = res     On return, a double which represents the final error value.

\author Mathieu Renouf.
*/


void lcp_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{

  int n;
  integer incx, incy;
  int i, iter;
  int itermax, ispeak;

  double err, a1, b1 , qs;

  double alpha, beta, rp, pMp;
  double den, num, tol;

  char NOTRANS = 'N';

  int *status;
  double *zz , *pp , *rr, *ww, *Mp;

  *info = 1;
  n     = *nn;
  incx  = 1;

  /*input*/

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  tol = dparamLCP[0];

  /*output*/

  iparamLCP[2] = 0;
  dparamLCP[1] = 0.0;

  qs = dnrm2_((integer *)&n , q , &incx);

  /*printf( " Norm: %g \n", qs );*/

  den = 1.0 / qs;

  /* Allocations */

  status = (int*)malloc(n * sizeof(int));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  incx = 1;

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

  dcopy_((integer *)&n , q , &incx , rr , &incy);

  a1 = -1.;
  b1 = -1.;

  dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , rr , &incy);

  /* Initialization of gradients */
  /* rr -> p and rr -> w */

  dcopy_((integer *)&n , rr , &incx , ww , &incy);
  dcopy_((integer *)&n , rr , &incx , pp , &incy);

  iter = 0.0;
  err  = 1.0 ;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Compute initial pMp */

    incx = 1;
    incy = 1;

    dcopy_((integer *)&n , pp , &incx , Mp , &incy);

    a1 = 1.0;
    b1 = 0.0;

    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , Mp , &incx , &b1 , w , &incy);

    pMp = ddot_((integer *)&n , pp , &incx , w  , &incy);

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
      free(status);

      iparamLCP[2] = iter;
      dparamLCP[1] = err;
      *info = 3;
      return;
    }

    rp  = ddot_((integer *)&n , pp , &incx , rr , &incy);

    alpha = rp / pMp;

    /*
     * Iterate prediction:
     * z' = z + alpha*p
     *
     */

    daxpy_((integer *)&n , &alpha , pp , &incx , z , &incy);

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

    dcopy_((integer *)&n , rr , &incx , w  , &incy);
    dcopy_((integer *)&n , q  , &incx , rr , &incy);

    a1 = -1.;
    b1 = -1.;

    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , rr , &incy);

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

    rp = ddot_((integer *)&n , ww , &incx, w , &incy);

    beta = -rp / pMp;

    dcopy_((integer *)&n , ww , &incx , pp , &incy);
    daxpy_((integer *)&n, &beta , zz , &incx , pp , &incy);

    /* **** Criterium convergence **** */

    qs   = -1.0;
    daxpy_((integer *)&n , &qs , rr , &incx , w , &incy);
    num = dnrm2_((integer *)&n, w , &incx);
    err = num * den;

  }

  iparamLCP[2] = iter;
  dparamLCP[1] = err;

  dcopy_((integer *)&n , rr , &incx , w , &incy);

  qs   = -1.0;
  dscal_((integer *)&n , &qs , w , &incx);


  *info = 1;
  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
    }
    else
    {
      printf(" Convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else if (err <= tol) *info = 0;

  free(Mp);
  free(status);

  free(ww);
  free(rr);
  free(pp);
  free(zz);

}
