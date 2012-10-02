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

#include "LA.h"
#include "NumericsOptions.h" // for global options
#include "FrictionContactProblem.h"
#include "SolverOptions.h"
#include "FrictionContact3D_compute_error.h"
#include "FrictionContact3D_projection.h"
#include "projectionOnCone.h"
#include "projectionOnCylinder.h"

#include <math.h>
#include <assert.h>
void FrictionContact3D_unitary_compute_and_add_error(double *z , double *w, double mu, double * error)
{

  double normUT;
  double rho = 1.0;

  double worktmp[3];
  /* Compute the modified local velocity */
  normUT = hypot(w[1], w[2]); // i.e sqrt(w[ic3p1]*w[ic3p1]+w[ic3p2]*w[ic3p2]);
  worktmp[0] = z[0] - rho * (w[0] + mu * normUT);
  worktmp[1] = z[1] - rho * w[1] ;
  worktmp[2] = z[2] - rho * w[2] ;
  projectionOnCone(worktmp, mu);
  worktmp[0] = z[0] -  worktmp[0];
  worktmp[1] = z[1] -  worktmp[1];
  worktmp[2] = z[2] -  worktmp[2];
  *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];

}
int FrictionContact3D_compute_error(
  FrictionContactProblem* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options, double * error)
{

  assert(problem);
  assert(z);
  assert(w);
  assert(error);

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double *mu = problem->mu;

  DCOPY(n , problem->q , incx , w , incy); // w <-q
  // Compute the current velocity
  prodNumericsMatrix(n, n, 1.0, problem->M, z, 1.0, w);

  *error = 0.;


  int ic, ic3;

  for (ic = 0, ic3 = 0 ; ic < nc ; ic++, ic3 += 3)
  {

    FrictionContact3D_unitary_compute_and_add_error(z + ic3, w + ic3, mu[ic], error);
  }
  *error = sqrt(*error);

  /* Computes error */
  double normq = DNRM2(n , problem->q , incx);
  *error = *error / (normq + 1.0);
  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - FrictionContact3D_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
    return 0;
}



int FrictionContact3D_compute_error_velocity(FrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions *options, double * error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("FrictionContact3D_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double *mu = problem->mu;
  double worktmp[3] = {0.0, 0.0, 0.0};
  double invmu = 0.0;
  DCOPY(n , problem->q , incx , z , incy); // z <-q

  // Compute the current reaction
  prodNumericsMatrix(n, n, 1.0, problem->M, w, 1.0, z);

  *error = 0.;
  double normUT = 0.0;
  double rho = 1.0;
  for (int ic = 0 ; ic < nc ; ic++)
  {
    /* Compute the modified local velocity */
    normUT = sqrt(w[ic * 3 + 1] * w[ic * 3 + 1] + w[ic * 3 + 2] * w[ic * 3 + 2]);
    worktmp[0] = w[ic * 3] + mu[ic] * normUT - rho * (z[ic * 3]);
    worktmp[1] = w[ic * 3 + 1] - rho * z[ic * 3 + 1] ;
    worktmp[2] = w[ic * 3 + 2] - rho * z[ic * 3 + 2] ;
    invmu = 1.0 / mu[ic];
    projectionOnCone(worktmp, invmu);
    normUT = sqrt(worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
    worktmp[0] = w[ic * 3] - (worktmp[0] - mu[ic] * normUT);
    worktmp[1] = w[ic * 3 + 1] -  worktmp[1];
    worktmp[2] = w[ic * 3 + 2] -  worktmp[2];
    *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }
  *error = sqrt(*error);

  /* Computes error */
  double normq = DNRM2(n , problem->q , incx);
  *error = *error / (normq + 1.0);
  if (*error > tolerance)
  {
    /*      if (verbose > 0) printf(" Numerics - FrictionContact3D_compute_error_velocity failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}


int FrictionContact3D_Tresca_compute_error(FrictionContactProblem* problem, double *z, double * w, double tolerance, SolverOptions * options, double* error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("FrictionContact3D_Tresca_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double worktmp[3];
  double R;
  DCOPY(n , problem->q , incx , w , incy); // w <-q
  // Compute the current velocity
  prodNumericsMatrix(n, n, 1.0, problem->M, z, 1.0, w);

  *error = 0.;
  double normUT;
  double rho = 1.0;
  for (int ic = 0 ; ic < nc ; ic++)
  {
    /* Compute the modified local velocity */
    normUT = sqrt(w[ic * 3 + 1] * w[ic * 3 + 1] + w[ic * 3 + 2] * w[ic * 3 + 2]);
    worktmp[0] = z[ic * 3] - rho * (w[ic * 3]);
    worktmp[1] = z[ic * 3 + 1] - rho * w[ic * 3 + 1] ;
    worktmp[2] = z[ic * 3 + 2] - rho * w[ic * 3 + 2] ;

    R = (options->dWork[ic]);
    projectionOnCylinder(worktmp, R);
    worktmp[0] = z[ic * 3] -  worktmp[0];
    worktmp[1] = z[ic * 3 + 1] -  worktmp[1];
    worktmp[2] = z[ic * 3 + 2] -  worktmp[2];
    *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }
  *error = sqrt(*error);

  /* Computes error */
  double normq = DNRM2(n , problem->q , incx);
  *error = *error / (normq + 1.0);
  if (*error > tolerance)
  {
    /*      if (verbose > 0) printf(" Numerics - FrictionContact3D_Tresca_compute_error failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}
