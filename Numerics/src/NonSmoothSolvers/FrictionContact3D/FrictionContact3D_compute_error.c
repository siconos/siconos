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
#include "FrictionContact_Problem.h"
#include "FrictionContact3D_projection.h"
#include <math.h>

int FrictionContact3D_compute_error(FrictionContact_Problem* problem, double *z , double *w, double tolerance, double * error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("FrictionContact3D_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double *mu = problem->mu;
  double worktmp[3];
  DCOPY(n , problem->q , incx , w , incy); // w <-q
  prod(n, n, 1.0, problem->M, z, 1.0, w);

  *error = 0.;
  double normUT;
  double rho = 1.0;
  for (int ic = 0 ; ic < nc ; ic++)
  {
    /* Compute the modified local velocity */
    normUT = sqrt(w[ic * 3 + 1] * w[ic * 3 + 1] + w[ic * 3 + 2] * w[ic * 3 + 2]);
    worktmp[0] = z[ic * 3] - rho * (w[ic * 3] + mu[ic] * normUT);
    worktmp[1] = z[ic * 3 + 1] - rho * w[ic * 3 + 1] ;
    worktmp[2] = z[ic * 3 + 2] - rho * w[ic * 3 + 2] ;
    projectionOnCone(worktmp, mu[ic]);
    worktmp[0] = z[ic * 3] -  worktmp[0];
    worktmp[1] = z[ic * 3 + 1] -  worktmp[1];
    worktmp[2] = z[ic * 3 + 2] -  worktmp[2];
    *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }
  *error = sqrt(*error);

  /* Computes error */
  double normq = DNRM2(n , problem->q , incx);
  *error = *error / normq;
  if (*error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - FrictionContact3D_compute_error failed: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}
