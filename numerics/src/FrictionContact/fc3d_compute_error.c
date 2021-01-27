/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/


#include "fc3d_compute_error.h"
#include <assert.h>                  // for assert
#include <float.h>                   // for DBL_EPSILON
#include <math.h>                    // for sqrt, fabs
#include <stddef.h>                  // for NULL
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "NumericsMatrix.h"          // for NM_prod_mv_3x3, NM_gemv
#include "SolverOptions.h"           // for SolverOptions
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"                   // for DEBUG_PRINTF, DEBUG_EXPR, DEBUG_...
#include "numerics_verbose.h"        // for numerics_error
#include "projectionOnCone.h"        // for projectionOnCone
#include "projectionOnCylinder.h"    // for projectionOnCylinder
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2


void fc3d_unitary_compute_and_add_error(double* restrict r, double* restrict u, double mu, double* restrict error, double * worktmp)
{

  //double normUT;
  //double worktmp[3];
  /* Compute the modified local velocity */
  /* worktmp[0] = r[0] - u[0] - mu *  hypot(u[1], u[2]); */
  worktmp[0] = r[0] -  u[0] - mu  * sqrt(u[1] * u[1] + u[2] * u[2]);
  worktmp[1] = r[1] -  u[1] ;
  worktmp[2] = r[2] -  u[2] ;
  projectionOnCone(worktmp, mu);
  worktmp[0] = r[0] -  worktmp[0];
  worktmp[1] = r[1] -  worktmp[1];
  worktmp[2] = r[2] -  worktmp[2];
  *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
}

int fc3d_compute_error(
  FrictionContactProblem* problem,
  double *z, double *w, double tolerance,
  SolverOptions * options, double norm, double * error)
{
  DEBUG_BEGIN("fc3d_compute_error(...)\n");
  assert(problem);
  assert(z);
  assert(w);
  assert(error);

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double *mu = problem->mu;

  /* Compute the current velocity */
  cblas_dcopy(n, problem->q, incx, w, incy);     // w <-q
  NM_prod_mv_3x3(n, n, problem->M, z, w); // w = Mz +q

  DEBUG_PRINTF("norm of the reaction %e\n", cblas_dnrm2(n, z, 1));
  DEBUG_PRINTF("norm of the velocity %e\n", cblas_dnrm2(n, w, 1));
  DEBUG_PRINTF("norm of q = %12.8e\n", norm);
  /* DEBUG_EXPR(NV_display(problem->q,n);); */
  /* DEBUG_EXPR(NV_display(w,n);); */
  /* DEBUG_EXPR(NV_display(z,n);); */

  *error = 0.;
  int ic, ic3;
  double worktmp[3];
  for(ic = 0, ic3 = 0 ; ic < nc ; ic++, ic3 += 3)
  {
    fc3d_unitary_compute_and_add_error(z + ic3, w + ic3, mu[ic], error, worktmp);
    /*DEBUG_PRINTF("absolute error = %12.8e contact =%i nc= %i\n", *error, ic, nc);*/
  }
  *error = sqrt(*error);
  DEBUG_PRINTF("absolute error in complementarity = %12.8e\n", *error);

  /* Compute relative error */
  double norm_r =cblas_dnrm2(n, z, 1);
  double norm_u =cblas_dnrm2(n, w, 1);
  double relative_scaling = fmax(norm, fmax(norm_r,norm_u)); 
  /* double relative_scaling = fmax(norm_r,norm_w); */
  /* double relative_scaling = norm; */

  if(fabs(relative_scaling) > DBL_EPSILON)
    *error /= relative_scaling;

  DEBUG_PRINTF("relative error in complementarity = %12.8e\n", *error);
  DEBUG_END("fc3d_compute_error(...)\n");
  if(*error > tolerance)
    return 1;

  return 0;
}



int fc3d_compute_error_velocity(FrictionContactProblem* problem, double *z, double *w, double tolerance,
                                SolverOptions *options, double * error)
{
  /* Checks inputs */
  if(problem == NULL || z == NULL || w == NULL)
    numerics_error("fc3d_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double *mu = problem->mu;
  double worktmp[3] = {0.0, 0.0, 0.0};
  double invmu = 0.0;
  cblas_dcopy(n, problem->q, incx, z, incy);     // z <-q

  // Compute the current reaction
  NM_gemv(1.0, problem->M, w, 1.0, z);

  *error = 0.;
  double normUT = 0.0;
  double rho = 1.0;
  for(int ic = 0 ; ic < nc ; ic++)
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
  double norm_q = cblas_dnrm2(n, problem->q, incx);
  *error = *error / (norm_q + 1.0);
  if(*error > tolerance)
  {
    /*      if (verbose > 0) printf(" Numerics - fc3d_compute_error_velocity failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}

void fc3d_Tresca_unitary_compute_and_add_error(double *z, double *w, double R, double * error, double * worktmp)
{
  /* Compute the modified local velocity */
  worktmp[0] = z[0] -  w[0];
  worktmp[1] = z[1] -  w[1] ;
  worktmp[2] = z[2] -  w[2] ;
  projectionOnCylinder(worktmp, R);
  worktmp[0] = z[0] -  worktmp[0];
  worktmp[1] = z[1] -  worktmp[1];
  worktmp[2] = z[2] -  worktmp[2];
  *error +=  worktmp[0] * worktmp[0]
             + worktmp[1] * worktmp[1]
             + worktmp[2] * worktmp[2];

}
int fc3d_Tresca_compute_error(FrictionContactProblem* problem,
                              double *z, double * w,
                              double tolerance, SolverOptions * options,
                              double norm,
                              double* error)
{
  /* Checks inputs */
  if(problem == NULL || z == NULL || w == NULL)
    numerics_error("fc3d_Tresca_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double R;
  cblas_dcopy(n, problem->q, incx, w, incy);     // w <-q
  // Compute the current velocity
  /* NM_gemv(1.0, problem->M, z, 1.0, w); */
  NM_prod_mv_3x3(n, n, problem->M, z, w);
  *error = 0.;
  int ic, ic3;
  double worktmp[3];
  for(ic = 0, ic3 = 0 ; ic < nc ; ic++, ic3 += 3)
  {
    R = (options->dWork[ic]);
    fc3d_Tresca_unitary_compute_and_add_error(z+ic3, w+ic3, R, error, worktmp);
  }
  *error = sqrt(*error);

  /* Computes error */
  DEBUG_PRINTF("norm = %12.8e\n", norm);
  if(fabs(norm) > DBL_EPSILON)
    *error /= norm;

  if(*error > tolerance)
  {
    /* if (verbose > 0) printf(" Numerics - fc3d_Tresca_compute_error failed: error = %g > tolerance = %g.\n",*error, tolerance);  */
    return 1;
  }
  else
    return 0;
}
