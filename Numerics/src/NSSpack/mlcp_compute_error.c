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
/*!\file mlcp_compute_error.c
 *
 * This function checks the validity of the vector z as a solution \n
 * of the MLCP : \n
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *   A u + Cv +a =0\\
 *   D u + Bv +b = w
 *   0 \le v \perp  w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
 * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
 * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
 * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
 * \author Vincent Acary form the routine  filter_result_MLCP.c of Pascal Denoyelle
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LA.h"
#include <math.h>

int mlcp_compute_error(int* nn, int* mm,  double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v,  int verbose, double *w, double *err)
{
  double error, normb;
  double errore, norma;
  double a1, b1;
  double *we;
  int i, incx, incy, n, m;
  int param = 1;
  n = *nn;
  m = *mm;

  incx = 1;
  incy = 1;

  a1 = 1.;
  b1 = 1.;
  we   = (double*)calloc(n, sizeof(double));
  DCOPY(n , a , incx , we , incy);  //  we <-- a
  DCOPY(m , b , incx , w , incy);  //  w <-- b


  // following int param, we recompute the product w = Du+BV +b and we = Au+CV +a
  // The test is then more severe if we compute w because it checks that the linear equation is satisfied

  if (param == 1)
  {
    DGEMV(LA_NOTRANS , m, n , a1 , D , m , u ,
          incx , b1 , w , incy);  // w <-- D*u+ w
    DGEMV(LA_NOTRANS , m , m , a1 , B , m , v ,
          incx , b1 , w , incy);  // w <-- B*v + w
  }



  DGEMV(LA_NOTRANS , n, n , a1 , A , n , u ,
        incx , b1 , we , incy);  // we <-- A*u+ we


  DGEMV(LA_NOTRANS , n , m , a1 , C , n , v ,
        incx , b1 , we , incy);  // we <-- C*v + we



  errore = 0.;
  errore =  DNRM2(n , we , incx);;

  error = 0.0;
  for (i = 0 ; i < m; i++)
  {

    if (v[i] < 0.0)
    {
      error += -v[i];
      if (w[i] < 0.0) error += v[i] * w[i];
    }
    if (w[i] < 0.0) error += -w[i];
    if ((v[i] > 0.0) && (w[i] > 0.0)) error += v[i] * w[i];
  }


  incx  = 1;
  normb = DNRM2(m , b , incx);
  norma = DNRM2(n , a , incx);





  if (error / normb >= errore / norma)
  {
    *err = error / (1.0 + normb);
  }
  else
  {
    *err = errore / (1.0 + norma);
  }
  free(we);
  if (verbose > 0) printf("Siconos/Numerics: mlcp_compute_error: Error evaluation = %g \n", *err);
  return 0;

}
