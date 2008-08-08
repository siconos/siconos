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
#include <math.h>


int FrictionContact2D_compute_error(FrictionContact_Problem* problem, double *z , double *w, double tolerance, double * error)
{

  int nc = problem->numberOfContacts;

  int n = nc * 2;

  int ic, iN, iT;

  double *mu = problem->mu;

  double tmp[3];

  double normT;

  /* Checks inputs */
  if (! problem || ! z || ! w)
    numericsError("FrictionContact2D_compute_error", "null input for problem and/or z and/or w");

  /* w = Mz + q computed in frictionContact2DLocalSolve */

  *error = 0.;

  /* Num. Methods For Nonsmooth Dynamics, A.13 P 528 */
  /* k* -) x _|_ y  (- k  <=>  x = projk(x-rho.y) */

  for (ic = 0, iN = 0, iT = 1 ; ic < nc ; ++ic , ++iN, ++iN, ++iT, ++iT)
  {
    /* Compute the modified local velocity */
    tmp[0] = z[iN] - (w[iN] + mu[ic] * fabs(w[iT]));  /* rho=1 */
    tmp[1] = z[iT] - w[iT];                     /* rho=1 */

    /* projection */
    normT = fabs(tmp[1]);
    if (mu[ic]*normT <= -tmp[0])
    {
      tmp[0] = 0.;
      tmp[1] = 0.;
    }
    else if (normT > mu[ic]*tmp[0])
    {
      /* solve([sqrt((r1-mu*ra)^2+(r0-ra)^2)=abs(mu*r0-r1)/sqrt(mu*mu+1)],[ra]) */
      tmp[0] = (mu[ic] * normT + tmp[0]) / (mu[ic] * mu[ic] + 1);
      tmp[1] = mu[ic] * tmp[0];
    }

    tmp[0] = z[iN] -  tmp[0];
    tmp[1] = z[iT] -  tmp[1];
    *error += tmp[0] * tmp[0] + tmp[1] * tmp[1];

  }
  *error = sqrt(*error);

  /* Computes error */
  *error /= DNRM2(n, problem->q, 1);

  if (*error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - FrictionContact2D_compute_error failed: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}
