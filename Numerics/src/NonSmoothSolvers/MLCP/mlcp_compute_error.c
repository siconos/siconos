/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
#include "LA.h"
#include "MLCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int mlcp_compute_error(MixedLinearComplementarity_Problem* problem, double *z, double *w, double tolerance, double * error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("mlcp_compute_error", "null input for problem and/or z and/or w");

  int param = 1;
  int n = problem->n; /* Equalities */
  int m = problem->m; /* Inequalities */
  int incx = 1, incy = 1;
  int size = n + m;

  /* Computation of w: depends on the way the problem is written */

  /* Problem in the form (M,q) */
  if (problem->problemType == 0)
  {
    if (problem->M == NULL)
      numericsError("mlcp_compute_error", "null input for M");

    /* Computes w = Mz + q */
    DCOPY(size , problem->q , incx , w , incy);
    prod(size, size, 1.0, problem->M, z, 1.0, w);

  }
  /* Problem in the form ABCD */
  else // if(problem->problemType == 1)
  {

    /* Checks inputs */
    if (problem->A == NULL || problem->B == NULL || problem->C == NULL  || problem->D == NULL);
    numericsError("mlcp_compute_error", "null input for A, B,C or D");

    /* Links to problem data */
    double *a = &problem->q[0];
    double *b = &problem->q[n];
    double *A = problem->A;
    double *B = problem->B;
    double *C = problem->C;
    double *D = problem->D;

    /* Compute "equalities" part, we = Au + Cv + a - Must be equal to 0 */
    DCOPY(n , a , incx , w , incy);  //  we = w[0..n-1] <-- a
    DGEMV(LA_NOTRANS , n, n , 1.0 , A , n , &z[0] , incx , 1.0 , w , incy);   // we <-- A*u + we
    DGEMV(LA_NOTRANS , n, m , 1.0 , C , n , &z[n] , incx , 1.0 , w , incy);   // we <-- C*v + we

    /* Computes part which corresponds to complementarity */
    double * wi = &w[n]; // No copy!!
    DCOPY(m , b , incx , wi , incy);  //  wi = w[n..m] <-- b
    // following int param, we recompute the product wi = Du+BV +b and we = Au+CV +a
    // The test is then more severe if we compute w because it checks that the linear equation is satisfied
    if (param == 1)
    {
      DGEMV(LA_NOTRANS , m, n , 1.0 , D , m , &z[0] , incx , 1.0 , wi , incy);   // wi <-- D*u+ w
      DGEMV(LA_NOTRANS , m , m , 1.0 , B , m , &z[n] , incx , 1.0 , wi , incy);  // wi <-- B*v + w
    }
  }

  /* Error on equalities part */
  double error_e = DNRM2(n , w , incx);;

  /* Checks complementarity (only for rows number n to size) */
  double error_i = 0.;
  double zi, wi;
  for (int i = n ; i < size ; i++)
  {
    zi = z[i];
    wi = w[i];
    if (zi < 0.0)
    {
      error_i += -zi;
      if (wi < 0.0) error_i += zi * wi;
    }
    if (wi < 0.0) error_i += -wi;
    if ((zi > 0.0) && (wi > 0.0)) error_i += zi * wi;
  }

  /* Computes error */
  double *q = problem->q;
  double normb = DNRM2(m , &q[n] , incx);
  double norma = DNRM2(n , &q[0] , incx);

  if (error_i / normb >= error_e / norma)
  {
    *error = error_i / (1.0 + normb);
  }
  else
  {
    *error = error_e / (1.0 + norma);
  }

  if (*error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - mlcp_compute_error failed: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
  {
    if (verbose > 0) printf("Siconos/Numerics: mlcp_compute_error: Error evaluation = %g \n", *error);
    return 0;
  }
}
