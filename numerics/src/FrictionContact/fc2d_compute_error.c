/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "fc2d_compute_error.h"
#include <float.h>                   // for DBL_EPSILON
#include <math.h>                    // for fabs, sqrt
#include <stdio.h>                   // for printf
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "NumericsMatrix.h"          // for NM_gemv
#include "numerics_verbose.h"        // for numerics_error, verbose
#include "SiconosBlas.h"                   // for cblas_dcopy

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"


#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

void fc2d_unitary_compute_and_add_error(double* restrict r, double* restrict u, double mu, double* restrict error, double * worktmp)
{
  /* Compute the modified local velocity */
  worktmp[0] = r[0] -  u[0] - mu  * fabs(u[1]);
  worktmp[1] = r[1] -  u[1] ;
  /* projection */
  double normT = fabs(worktmp[1]);
  if(mu*normT <= -worktmp[0])
  {
    worktmp[0] = 0.;
    worktmp[1] = 0.;
  }
  else if(normT > mu*worktmp[0])
  {
    worktmp[0] = (mu * normT + worktmp[0]) / (mu * mu + 1.0);
    worktmp[1] = mu * worktmp[0] * SGN(worktmp[1]);
  }
  
  worktmp[0] = r[0] -  worktmp[0];
  worktmp[1] = r[1] -  worktmp[1];
  *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1];
}



/* static int fc2d_compute_error_old( */
/*   FrictionContactProblem* problem, */
/*   double *z, */
/*   double *w, */
/*   double tolerance, */
/*   double norm, */
/*   double * error) */
/* { */
/*   DEBUG_BEGIN("fc2d_compute_error(...)\n"); */
/*   /\* Checks inputs *\/ */
/*   if(! problem || ! z || ! w) */
/*     numerics_error("fc2d_compute_error", "null input for problem and/or z and/or w"); */

/*   int nc = problem->numberOfContacts; */

/*   int n = nc * 2; */

/*   int ic, iN; */

/*   double *mu = problem->mu; */

/*   double tmp[2]; */

/*   cblas_dcopy(n, problem->q, 1, w, 1); // w <-q */
/*   NM_gemv(1.0, problem->M, z, 1.0, w); */

/*   *error = 0.; */

/*   for(ic = 0, iN = 0 ; ic < nc ; ++ic, ++iN, ++iN) */
/*   { */
/*     fc2d_unitary_compute_and_add_error( &z[iN], &w[iN], mu[ic], error, tmp); */
/*   } */

/*   DEBUG_PRINTF("error=%e\n", *error ); */

/*   *error = sqrt(*error); */
/*   if(fabs(norm) > DBL_EPSILON) */
/*     *error /= norm; */

/*   DEBUG_END("fc2d_compute_error(...)\n"); */
/*   if(*error > tolerance) */
/*   { */
/*     numerics_printf_verbose(1, "-- FC2D - warning. fc2d_compute_error = %g > %g (= tolerance)", *error, tolerance); */
/*     return 1; */
/*   } */
/*   else */
/*     return 0; */
/* } */

int fc2d_compute_error(
  FrictionContactProblem* problem,
  double *z,
  double *w,
  double tolerance,
  double norm,
  double * error)
{
  DEBUG_BEGIN("fc2d_compute_error(...)\n");
  /* Checks inputs */
  assert(z);
  assert(w);
  assert(problem);
  assert(error);

  int nc = problem->numberOfContacts;
  int n = nc * 2;

  double *mu = problem->mu;
  double tmp[2];

  cblas_dcopy(n, problem->q, 1, w, 1); // w <-q
  NM_gemv(1.0, problem->M, z, 1.0, w); // w = W z +q is assumed to be compute at machine accuracy.

  *error = 0.;
  for(int ic = 0, iN = 0 ; ic < nc ; ++ic, ++iN, ++iN)
    fc2d_unitary_compute_and_add_error( &z[iN], &w[iN], mu[ic], error, tmp);
  DEBUG_PRINTF("error=%e\n", *error );

  *error = sqrt(*error);

  double norm_r = cblas_dnrm2(n, z, 1);
  double norm_w = cblas_dnrm2(n, w, 1);
  double relative_scaling = fmax(norm_r, fmax(norm_w,norm));
  if(fabs(relative_scaling) > DBL_EPSILON)
    *error /= relative_scaling;

  DEBUG_END("fc2d_compute_error(...)\n");
  if(*error > tolerance)
  {
    numerics_printf_verbose(1, "-- FC2D - warning. fc2d_compute_error = %g > %g (= tolerance)", *error, tolerance);
    return 1;
  }
  else
    return 0;
}
