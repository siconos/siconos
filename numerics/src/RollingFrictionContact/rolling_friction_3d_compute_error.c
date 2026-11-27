// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "rolling_friction_3d_compute_error.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs

#include "NumericsMatrix.h"                 // for NM_gemv
#include "RollingFrictionContactProblem.h"  // for RollingFrictionContactPro...
#include "SiconosBlas.h"                    // for cblas_dcopy
#include "convex_set.h"                     // for ConvexSet, convex_set_natural_map_block
#include "numerics_errors.h"

// #define DEBUG_MESSAGES                                       //
#include "siconos_debug.h"  // for DEBUG_EXPR, DEBUG_PRINTF
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

void rolling_friction_3d_unitary_compute_and_add_error(double r[5], double u[5], double mu,
                                                       double mur, double* restrict error,
                                                       double* worktmp) {
  /* natural map r - P_K(r - u_mod), u_mod = u + (mu ||u_T|| + mur ||u_R||) e_N the modified
     velocity, K the rolling friction cone */
  worktmp[0] =
      u[0] + mu * sqrt(u[1] * u[1] + u[2] * u[2]) + mur * sqrt(u[3] * u[3] + u[4] * u[4]);
  for (int k = 1; k < 5; k++) worktmp[k] = u[k];
  ConvexSet K = convex_set(CONVEX_SET_ROLLING_CONE, 5, 1, &mu, &mur);
  *error += convex_set_natural_map_block(&K, 0, r, worktmp, NULL);
}

int rolling_friction_3d_compute_error(RollingFrictionContactProblem* problem, double* reaction,
                                      double* velocity, double tolerance,
                                      SolverOptions* options, double norm, double* error) {
  DEBUG_BEGIN("rolling_friction_3d_compute_error(...)\n");
  CHECK_NULL(problem);
  CHECK_NULL(reaction);
  CHECK_NULL(velocity);
  CHECK_NULL(error);

  /* Computes velocity = Mreaction + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 5;
  double* mu = problem->mu;
  double* mur = problem->mu_r;

  cblas_dcopy(n, problem->q, incx, velocity, incy);  // velocity <-q
  // Compute the current velocity
  NM_gemv(1.0, problem->M, reaction, 1.0, velocity);
  DEBUG_EXPR(NV_display(problem->q, n););
  DEBUG_EXPR(NV_display(velocity, n););
  DEBUG_EXPR(NV_display(reaction, n););

  *error = 0.;
  int ic, ic5;
  double worktmp[5];
  for (ic = 0, ic5 = 0; ic < nc; ic++, ic5 += 5) {
    rolling_friction_3d_unitary_compute_and_add_error(reaction + ic5, velocity + ic5, mu[ic],
                                                      mur[ic], error, worktmp);
    DEBUG_PRINTF("squared absolute error = %12.8e contact =%i nc= %i\n", *error, ic, nc);
  }
  *error = sqrt(*error);
  DEBUG_PRINTF("absolute error = %12.8e\n", *error);
  /* Compute relative error with respect to norm */
  DEBUG_PRINTF("norm = %12.8e\n", norm);

  norm = fmax(norm, cblas_dnrm2(n, reaction, 1));
  norm = fmax(norm, cblas_dnrm2(n, velocity, 1));
  if (fabs(norm) > DBL_EPSILON) *error /= norm;
  /* *error = *error / (norm + 1.0); old version */
  DEBUG_PRINTF("relative error = %12.8e\n", *error);

  DEBUG_END("rolling_friction_3d_compute_error(...)\n");
  if (*error > tolerance) return 1;

  return 0;
}
