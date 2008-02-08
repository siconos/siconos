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

#include "LA.h"
#include "Numerics_Options.h" // for global options
#include "LCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int filter_result_LCP(LinearComplementarity_Problem* problem, double *z , double *w, double tolerance)
{

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("filter_result_LCP", "null input for problem and/or z and/or w");

  if (problem->M->storageType == 1)
    numericsError("filter_result_LCP", "Not yet implemented for sparse storage");

  double error, normq;
  int i, incx, incy;
  double zi, wi;

  /* get var. from problem (pointer links, no copy!)*/
  double *q = problem->q;
  double *M = problem->M->matrix0;
  int n = problem->size;

  incx = 1;
  incy = 1;
  DCOPY(n , q , incx , w , incy);
  DGEMV(LA_NOTRANS , n , n , 1.0 , M , n , z , incx , 1.0 , w , incy);

  error = 0.;
  for (i = 0 ; i < n ; i++)
  {
    zi = z[i];
    wi = w[i];
    if (zi < 0.0)
    {
      error += -zi;
      if (wi < 0.0) error += zi * wi;
    }
    if (wi < 0.0) error += -wi;
    if ((zi > 0.0) && (wi > 0.0)) error += zi * wi;
  }

  incx  = 1;
  normq = DNRM2(n , q , incx);

  error = error / normq;
  if (error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - filter_result_LCP failed: error = %g > tolerance = %g.\n", error, tolerance);
    return 1;
  }
  else
    return 0;
}
