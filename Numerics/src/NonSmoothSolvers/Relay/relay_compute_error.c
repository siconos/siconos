/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "LA.h"
#include "NumericsOptions.h" // for global options
#include "RelayProblem.h"


int projectiononbox(double *z , double *lb, double * ub, int n)
{

  for (int i = 0 ; i < n ; i++)
  {
    if (z[i] < lb[i]) z[i] = lb[i];
    else if (z[i] > ub[i]) z[i] = ub[i];
  }
  return 0;
}

int relay_compute_error(RelayProblem* problem, double *z , double *w, double tolerance, double * error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("relay_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int n = problem->size;
  DCOPY(n , problem->q , 1 , w , 1);  // w <-q
  prodNumericsMatrix(n, n, 1.0, problem->M, z, 1.0, w);
  double * ztmp = (double*)malloc(n * sizeof(double));
  DCOPY(n , z , 1 , ztmp, 1);  // ztmp <-z

  double rho = -1.0;
  DAXPY(n, rho, w, 1, ztmp, 1);    //ztmp <- ztmp - rho w

  projectiononbox(ztmp , problem->lb, problem->ub, n);
  DAXPY(n, -1, z, 1, ztmp, 1);    //ztmp <- ztmp -z


  *error = DNRM2(n , ztmp , 1);


  /* Computes error */
  double normq = DNRM2(n , problem->q , 1);
  *error = *error / (normq + 1.0);
  free(ztmp);
  if (*error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - relay_compute_error: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}
