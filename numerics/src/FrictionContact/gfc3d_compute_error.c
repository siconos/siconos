/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include "gfc3d_compute_error.h"

#include "NumericsOptions.h" // for global options
#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"
#include "projectionOnCone.h"
#include "SiconosLapack.h"
#include <math.h>
#include <assert.h>
#include <float.h>
#include "sanitizer.h"

int gfc3d_compute_error(GlobalFrictionContactProblem* problem, double* restrict reaction , double* restrict velocity, double* restrict globalVelocity, double tolerance, double* restrict error)
{

  /* Checks inputs */
  if (problem == NULL || reaction == NULL || velocity == NULL || globalVelocity == NULL)
    numericsError("gfc3d_compute_error", "null input");

  gfc3d_init_workspace(problem);
  NumericsMatrix* factorized_M = problem->workspace->factorized_M;
  double* globalVelocitytmp = problem->workspace->globalVelocity;

  /* Computes error = dnorm2( GlobalVelocity -M^-1( q + H reaction)*/
  int nc = problem->numberOfContacts;
  int m = nc * 3;
  int n = problem->M->size0;
  double *mu = problem->mu;
  double *q = problem->q;
  NumericsMatrix *H = problem->H;

  cblas_dcopy_msan(n, q, 1, globalVelocitytmp, 1);

  NM_gemv(1.0, H, reaction, 1.0, globalVelocitytmp);

  CHECK_RETURN(!NM_gesv_expert(factorized_M, globalVelocitytmp, true));

  cblas_daxpy(n , -1.0 , globalVelocity , 1 , globalVelocitytmp, 1);

  *error =   cblas_dnrm2(n , globalVelocitytmp , 1);

  cblas_dcopy(m, problem->b, 1, velocity, 1);
  NM_tgemv(1, H, globalVelocity, 1, velocity);

  double worktmp[3];
  double normUT;
  double rho = 1.0;
  for (int ic = 0 ; ic < nc ; ic++)
  {
    /* Compute the modified local velocity */
    normUT = sqrt(velocity[ic * 3 + 1] * velocity[ic * 3 + 1] + velocity[ic * 3 + 2] * velocity[ic * 3 + 2]);
    worktmp[0] = reaction[ic * 3] - rho * (velocity[ic * 3] + mu[ic] * normUT);
    worktmp[1] = reaction[ic * 3 + 1] - rho * velocity[ic * 3 + 1] ;
    worktmp[2] = reaction[ic * 3 + 2] - rho * velocity[ic * 3 + 2] ;
    projectionOnCone(worktmp, mu[ic]);
    worktmp[0] = reaction[ic * 3] -  worktmp[0];
    worktmp[1] = reaction[ic * 3 + 1] -  worktmp[1];
    worktmp[2] = reaction[ic * 3 + 2] -  worktmp[2];
    *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }
  *error = sqrt(*error);

  /* Computes error */
  double normq = cblas_dnrm2(n , problem->q , 1);
  *error = *error / (normq + 1.0);

  if (*error > tolerance)
  {
    /*       if (verbose > 0) printf(" Numerics - gfc3d_compute_error failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}
