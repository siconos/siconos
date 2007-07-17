/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
/*!\file lcp_compute_error.c
 *
 * This function checks the validity of the vector z as a solution \n
 * of the LCP : \n
 * \f$
 *    0 \le z \perp Mz + q \ge 0
 * \f$
 * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
 * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
 * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
 * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
 * \author Vincent Acary form the routine  filter_result_LCP.c of Pascal Denoyelle
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blaslapack.h"
#include <math.h>

int lcp_compute_error(int n, double *vec , double *q , double *z , int verbose, double *w, double *err)
{
  double error, normq;
  double a1, b1;
  int i, incx, incy;
  int param = 1;

  char NOTRANS = 'N';

  incx = 1;
  incy = 1;
  dcopy_((integer *)&n , q , (integer *)&incx , w , (integer *)&incy);

  a1 = 1.;
  b1 = 1.;


  // following int param, we recompute the product w = M*z+q
  // The test is then more severe if we compute w because it checks that the linear equation is satisfied

  if (param == 1)
  {
    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z ,
           (integer *)&incx , &b1 , w , (integer *)&incy);
  }



  error = 0.;
  for (i = 0 ; i < n ; i++)
  {
    if (z[i] < 0.0)
    {
      error += -z[i];
      if (w[i] < 0.0) error += z[i] * w[i];
    }
    if (w[i] < 0.0) error += -w[i];
    if ((z[i] > 0.0) && (w[i] > 0.0)) error += z[i] * w[i];
  }

  incx  = 1;
  normq = dnrm2_((integer *)&n , q , (integer *)&incx);

  *err = error / normq;

  if (verbose > 0) printf("Siconos/Numerics: lcp_compute_error: Error evaluation = %g \n", *err);
  return 0;

}
