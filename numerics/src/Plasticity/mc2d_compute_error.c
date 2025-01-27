/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "mc2d_compute_error.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs
#include <stddef.h>  // for NULL

#include "MohrCoulomb2DProblem.h"  // for Mohrcoulomb2dproblem
#include "NumericsMatrix.h"        // for NM_prod_mv_3x3, NM_gemv
#include "SolverOptions.h"         // for SolverOptions
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "numerics_verbose.h"      // for numerics_error
#include "projectionOnCone.h"      // for projectionOnCone
#include "projectionOnCylinder.h"  // for projectionOnCylinder
#include "siconos_debug.h"         // for DEBUG_PRINTF, DEBUG_EXPR, DEBUG_...
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif
#include "SiconosBlas.h"  // for cblas_dcopy, cblas_dnrm2

void mc2d_unitary_compute_and_add_error(double *restrict r, double *restrict u, double eta,
                                        double theta, double *restrict error, double *worktmp) {
  // double normUT;
  // double worktmp[3];
  /* Compute the modified local velocity */
  /* worktmp[0] = r[0] - u[0] - eta *  hypot(u[1], u[2]); */
  worktmp[0] = r[0] - u[0] - theta * sqrt(u[1] * u[1] + u[2] * u[2]);
  worktmp[1] = r[1] - u[1];
  worktmp[2] = r[2] - u[2];
  projectionOnCone(worktmp, eta);
  worktmp[0] = r[0] - worktmp[0];
  worktmp[1] = r[1] - worktmp[1];
  worktmp[2] = r[2] - worktmp[2];
  *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
}

int mc2d_compute_error(MohrCoulomb2DProblem *problem, double *z, double *w, double tolerance,
                       SolverOptions *options, double norm, double *error) {
  DEBUG_BEGIN("mc2d_compute_error(...)\n");
  assert(problem);
  assert(z);
  assert(w);
  assert(error);

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfCones;
  int n = nc * 3;
  double *eta = problem->eta;
  double *theta = problem->theta;

  /* Compute the current velocity */
  cblas_dcopy(n, problem->q, incx, w, incy);  // w <-q
  NM_prod_mv_3x3(n, n, problem->M, z, w);     // w = Mz +q

  DEBUG_PRINTF("norm of the reaction %e\n", cblas_dnrm2(n, z, 1));
  DEBUG_PRINTF("norm of the velocity %e\n", cblas_dnrm2(n, w, 1));
  DEBUG_PRINTF("norm of q = %12.8e\n", norm);
  /* DEBUG_EXPR(NV_display(problem->q,n);); */
  /* DEBUG_EXPR(NV_display(w,n);); */
  /* DEBUG_EXPR(NV_display(z,n);); */

  *error = 0.;
  int ic, ic3;
  double worktmp[3];
  for (ic = 0, ic3 = 0; ic < nc; ic++, ic3 += 3) {
    mc2d_unitary_compute_and_add_error(z + ic3, w + ic3, eta[ic], theta[ic], error, worktmp);
    /*DEBUG_PRINTF("absolute error = %12.8e contact =%i nc= %i\n", *error, ic, nc);*/
  }
  *error = sqrt(*error);
  DEBUG_PRINTF("absolute error in complementarity = %12.8e\n", *error);

  /* Compute relative error */
  double norm_r = cblas_dnrm2(n, z, 1);
  double norm_u = cblas_dnrm2(n, w, 1);
  double relative_scaling = fmax(norm, fmax(norm_r, norm_u));
  /* double relative_scaling = fmax(norm_r,norm_w); */
  /* double relative_scaling = norm; */

  //if (fabs(relative_scaling) > DBL_EPSILON) *error /= relative_scaling;

  DEBUG_PRINTF("relative error in complementarity = %12.8e\n", *error);
  DEBUG_END("mc2d_compute_error(...)\n");
  if (*error > tolerance) return 1;

  return 0;
}
