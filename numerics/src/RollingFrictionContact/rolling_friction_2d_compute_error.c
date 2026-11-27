// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "rolling_friction_2d_compute_error.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs

#include "NumericsMatrix.h"                 // for NM_gemv
#include "RollingFrictionContactProblem.h"  // for RollingFrictionContactPro...
#include "SiconosBlas.h"                    // for cblas_dcopy
#include "convex_set.h"                     // for ConvexSet, convex_set_natural_map_block
#include "numerics_errors.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"  // for DEBUG_EXPR, DEBUG_PRINTF
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

void rolling_friction_2d_unitary_compute_and_add_error(double r[3], double u[3], double mu,
                                                       double mur, double* restrict error,
                                                       double* worktmp) {
  /* natural map r - P_K(r - u_mod), u_mod = u + (mu |u_T| + mur |u_R|) e_N the modified
     velocity, K the 2D rolling friction cone */
  worktmp[0] = u[0] + mu * fabs(u[1]) + mur * fabs(u[2]);
  worktmp[1] = u[1];
  worktmp[2] = u[2];
  ConvexSet K = convex_set(CONVEX_SET_ROLLING_CONE, 3, 1, &mu, &mur);
  *error += convex_set_natural_map_block(&K, 0, r, worktmp, NULL);
}

int rolling_friction_2d_compute_error(RollingFrictionContactProblem* problem, double* reaction,
                                      double* velocity, double tolerance,
                                      SolverOptions* options, double norm, double* error) {
  DEBUG_BEGIN("rolling_friction_2d_compute_error(...)\n");
  CHECK_NULL(problem);
  CHECK_NULL(reaction);
  CHECK_NULL(velocity);
  CHECK_NULL(error);

  /* Computes velocity = Mreaction + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double* mu = problem->mu;
  double* mur = problem->mu_r;

  cblas_dcopy(n, problem->q, incx, velocity, incy);  // velocity <-q
  // Compute the current velocity
  NM_prod_mv_3x3(n, n, problem->M, reaction, velocity);
  /* NM_gemv(1.0, problem->M, reaction, */
  /*         1.0, */
  /*         velocity); */
  DEBUG_EXPR(NV_display(problem->q, n););
  DEBUG_EXPR(NV_display(velocity, n););
  DEBUG_EXPR(NV_display(reaction, n););

  *error = 0.;
  int ic, ic3;
  double worktmp[3];
  for (ic = 0, ic3 = 0; ic < nc; ic++, ic3 += 3) {
    rolling_friction_2d_unitary_compute_and_add_error(reaction + ic3, velocity + ic3, mu[ic],
                                                      mur[ic], error, worktmp);
    DEBUG_PRINTF("squared absolute error = %12.8e contact =%i nc= %i\n", *error, ic, nc);
  }
  *error = sqrt(*error);
  DEBUG_PRINTF("absolute error = %12.8e\n", *error);
  /* Compute relative error with respect to norm */
  DEBUG_PRINTF("norm = %12.8e\n", norm);
  if (fabs(norm) > DBL_EPSILON) *error /= norm;
  /* *error = *error / (norm + 1.0); old version */
  DEBUG_PRINTF("relative error = %12.8e\n", *error);
  DEBUG_END("rolling_friction_2d_compute_error(...)\n");
  if (*error > tolerance) return 1;

  return 0;
}
