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
#include "Numerics_Options.h" // for global options
#include "LinearComplementarity_Problem.h"

int lcp_compute_error(LinearComplementarity_Problem* problem, double *z , double *w, double tolerance, double * error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("lcp_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int n = problem->size;
  DCOPY(n , problem->q , incx , w , incy);
  prod(n, n, 1.0, problem->M, z, 1.0, w);

  /* Checks complementarity */
  *error = 0.;
  double zi, wi;
  for (int i = 0 ; i < n ; i++)
  {
    zi = z[i];
    wi = w[i];
    if (zi < 0.0)
    {
      *error += -zi;
      if (wi < 0.0) *error += zi * wi;
    }
    if (wi < 0.0) *error += -wi;
    if ((zi > 0.0) && (wi > 0.0)) *error += zi * wi;
  }

  /* Computes error */
  double normq = DNRM2(n , problem->q , incx);
  *error = *error / normq;
  if (*error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - lcp_compute_error failed: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}
