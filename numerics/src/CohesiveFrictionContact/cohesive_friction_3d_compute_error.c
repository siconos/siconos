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

#include "cohesive_friction_3d_compute_error.h"

#include <assert.h>
#include <float.h>
#include <math.h>

#include "NumericsMatrix.h"
#include "SiconosBlas.h"
#include "fc3d_compute_error.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "projectionOnDisk.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"

void cohesive_friction_3d_unitary_compute_and_add_error(double r[3], double u[3], double mu,
                                                        double* error, double worktmp[3]) {
  /* Compute the modified local velocity: w = r - u */
  worktmp[0] = r[0] - u[0];
  worktmp[1] = r[1] - u[1];
  worktmp[2] = r[2] - u[2];

  /* Project onto the friction cone */
  projectionOnDisk(&worktmp[1], mu);

  /* Compute residual: r - P_C(r - u) */
  worktmp[0] = 0.0;
  worktmp[1] = r[1] - worktmp[1];
  worktmp[2] = r[2] - worktmp[2];

  /* Accumulate squared error */
  *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
}

int cohesive_friction_3d_compute_error(CohesiveFrictionContactProblem* problem,
                                       double* reaction, double* velocity, double tolerance,
                                       SolverOptions* options, double norm, double* error) {
  DEBUG_BEGIN("cohesive_friction_3d_compute_error(...)\n");

  CHECK_NULL(problem);
  CHECK_NULL(reaction);
  CHECK_NULL(velocity);
  CHECK_NULL(error);
  CHECK_MATRIX(problem->M);
  CHECK_NULL(problem->q);
 
  int nc = problem->numberOfContacts;
  int ncoh = problem->numberOfCohesivePoints;

  int dim = problem->dimension;
  int n = dim * (nc + ncoh);
  double* mu = problem->mu;

  *error = 0.0;
  double worktmp[3];
  int incx = 1, incy = 1;

  /* Compute the current velocity */
  cblas_dcopy(n, problem->q, incx, velocity, incy);      // w <-q
  NM_prod_mv_3x3(n, n, problem->M, reaction, velocity);  // w = Mz +q

  /* Loop over all contacts and accumulate error */
  for (int ic = 0, ic3 = 0; ic < nc; ic++, ic3 += dim) {
    fc3d_unitary_compute_and_add_error(reaction + ic3, velocity + ic3, mu[ic], error, worktmp);
  }
  /* Loop over all cohesive points and accumulate error */
  for (int ic = nc, ic3 = 3 * nc; ic < nc + ncoh; ic++, ic3 += dim) {
    cohesive_friction_3d_unitary_compute_and_add_error(reaction + ic3, velocity + ic3,
                                                       problem->c_t[ic - nc], error, worktmp);
  }

  *error = sqrt(*error);
  DEBUG_PRINTF("absolute error = %12.8e\n", *error);

  /* Compute relative error with proper normalization */
  double norm_r = cblas_dnrm2(n, reaction, 1);
  double norm_u = cblas_dnrm2(n, velocity, 1);
  
  /* DEBUG_PRINTF("norm_r = %2.4e\n", norm_r); */
  /* DEBUG_PRINTF("norm_u = %2.4e\n", norm_u); */
    
  double relative_scaling = fmax(norm, fmax(norm_r, norm_u));

  if (fabs(relative_scaling) > DBL_EPSILON) {
    *error /= relative_scaling;
  }

  DEBUG_PRINTF("relative error = %12.8e\n", *error);
  DEBUG_END("cohesive_friction_3d_compute_error(...)\n");

  return (*error > tolerance) ? 1 : 0;
}
