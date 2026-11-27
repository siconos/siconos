// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "fc3d_compute_error.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs
#include <stddef.h>  // for NULL
#include <stdlib.h>  // for malloc, free

#include "FrictionContactProblem.h"  // IWYU pragma: keep - For FrictionContactProblem
#include "FrictionContact_options.h"
#include "NumericsMatrix.h"     // for NM_prod_mv_3x3, NM_gemv
#include "SolverOptions.h"      // IWYU pragma: keep - for SolverOptions
#include "error_computation.h"  // for error_compute_alloc_pair
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "convex_set.h"         // for ConvexSet, convex_set_natural_map_block
#include "natural_map_error.h"  // for natural_map_error
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "projection_soc.h"  // for soc_dual_projection_3d
#include "siconos_debug.h"   // for DEBUG_PRINTF, DEBUG_EXPR, DEBUG_...
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif
#include "SiconosBlas.h"  // for cblas_dcopy, cblas_dnrm2

void fc3d_unitary_compute_and_add_error(const double r[3], const double u[3], double mu,
                                        double* restrict error, double worktmp[3]) {
  /* natural map r - P_K(r - u_mod), u_mod = u + mu ||u_T|| e_N the modified velocity */
  worktmp[0] = u[0] + mu * sqrt(u[1] * u[1] + u[2] * u[2]);
  worktmp[1] = u[1];
  worktmp[2] = u[2];
  ConvexSet K = convex_set(CONVEX_SET_SOC, 3, 1, &mu, NULL);
  *error += convex_set_natural_map_block(&K, 0, r, worktmp, NULL);
}

void fc3d_compute_modified_local_velocity(int nc, const double* u, const double* mu,
                                          double* u_modified, double* work) {
  int alloc_work = 0;
  if (work == NULL) {
    work = (double*)malloc(nc * sizeof(double));
    alloc_work = 1;
  }

  /* Tangential norms of the local velocity, all contacts at once.
     A plain loop is faster than per-block BLAS calls for these short vectors. */
  for (int i = 0; i < nc; ++i) {
    const double t1 = u[3 * i + 1];
    const double t2 = u[3 * i + 2];
    work[i] = sqrt(t1 * t1 + t2 * t2);
  }

  for (int i = 0; i < nc; ++i) {
    const int pos = 3 * i;
    u_modified[pos] = u[pos] + mu[i] * work[i];
    u_modified[pos + 1] = u[pos + 1];
    u_modified[pos + 2] = u[pos + 2];
  }

  if (alloc_work) free(work);
}

double fc3d_compute_complementarity_error(int nc, const double* r, const double* u,
                                          const double* mu, int norm, int scope, int reduction,
                                          double alpha, double* work, double* tmp) {
  /* modified local velocity u_mod in tmp, then the scaled natural-map error on the stack of
     friction cones (natural_map_error.h) */
  fc3d_compute_modified_local_velocity(nc, u, mu, tmp, work);
  NaturalMapErrorScaling s = {norm, scope, reduction, alpha};
  ConvexSet K = convex_set(CONVEX_SET_SOC, 3, nc, mu, NULL);
  return natural_map_error(&s, &K, r, tmp);
}

void fc3d_unitary_compute_dual_and_add_error(double* restrict r, double* restrict u, double mu,
                                             double* restrict error, double* worktmp) {
  worktmp[0] = u[0] + mu * sqrt(u[1] * u[1] + u[2] * u[2]) - r[0];
  worktmp[1] = u[1] - r[1];
  worktmp[2] = u[2] - r[2];
  soc_dual_projection_3d(worktmp, mu);
  worktmp[0] = u[0] + mu * sqrt(u[1] * u[1] + u[2] * u[2]) - worktmp[0];
  worktmp[1] = u[1] - worktmp[1];
  worktmp[2] = u[2] - worktmp[2];
  *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
}

int fc3d_compute_error(FrictionContactProblem* problem, double* z, double* w, double tolerance,
                       SolverOptions* options, double* error) {
  DEBUG_BEGIN("fc3d_compute_error(...)\n");
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(w);
  CHECK_NULL(error);
  CHECK_MATRIX(problem->M);
  CHECK_NULL(problem->q);
  CHECK_NULL(problem->mu);

  /* Computes w = Mz + q */
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double* mu = problem->mu;

  /* Compute the current velocity */
  cblas_dcopy(n, problem->q, 1, w, 1);     // w <- q
  NM_prod_mv_3x3(n, n, problem->M, z, w);  // w = Mz + q

  DEBUG_PRINTF("norm of the reaction %e\n", cblas_dnrm2(n, z, 1));
  DEBUG_PRINTF("norm of the velocity %e\n", cblas_dnrm2(n, w, 1));

  const NaturalMapErrorScaling defaults = {SICONOS_FRICTION_3D_SCALING_EUCLIDEAN,
                                           SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
                                           SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, 1.0};
  NaturalMapErrorScaling sc = natural_map_error_scaling(options, defaults);
  double* comp_work = NULL;
  double* comp_tmp = NULL;
  if (error_compute_alloc_pair(nc, n, &comp_work, &comp_tmp)) {
    return numerics_error("fc3d_compute_error", "bad alloc");
  }
  *error = fc3d_compute_complementarity_error(nc, z, w, mu, sc.norm, sc.scope, sc.reduction,
                                              sc.alpha, comp_work, comp_tmp);
  free(comp_work);
  free(comp_tmp);
  DEBUG_PRINTF("scaled error in complementarity = %12.8e\n", *error);

  DEBUG_END("fc3d_compute_error(...)\n");
  if (*error > tolerance) return 1;

  return 0;
}

int fc3d_compute_error_norm_infinity_conic(FrictionContactProblem* problem, double* z,
                                           double* w, double tolerance, SolverOptions* options,
                                           double norm, double* error, int on_dual_cone) {
  DEBUG_BEGIN("fc3d_compute_error_norm_infinity_conic(...)\n");
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(w);
  CHECK_NULL(error);

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double* mu = problem->mu;
  double error_unitary = 0.;

  /* Compute the current velocity */
  cblas_dcopy(n, problem->q, incx, w, incy);  // w <-q
  NM_prod_mv_3x3(n, n, problem->M, z, w);     // w = Mz +q

  *error = -1.;
  int ic, ic3;
  double worktmp[3];
  for (ic = 0, ic3 = 0; ic < nc; ic++, ic3 += 3) {
    error_unitary = 0.;
    if (on_dual_cone) {
      fc3d_unitary_compute_dual_and_add_error(z + ic3, w + ic3, mu[ic], &error_unitary,
                                              worktmp);
    } else {
      fc3d_unitary_compute_and_add_error(z + ic3, w + ic3, mu[ic], &error_unitary, worktmp);
    }

    *error = fmax(*error, error_unitary);
  }
  *error = sqrt(*error);

  return 0;
}

int fc3d_compute_error_velocity(FrictionContactProblem* problem, double* z, double* w,
                                double tolerance, SolverOptions* options, double* error) {
  (void)options; /* unused */
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    return numerics_error("fc3d_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double* mu = problem->mu;
  double worktmp[3] = {0.0, 0.0, 0.0};
  double invmu = 0.0;
  cblas_dcopy(n, problem->q, incx, z, incy);  // z <-q

  // Compute the current reaction
  NM_gemv(1.0, problem->M, w, 1.0, z);

  *error = 0.;
  double normUT = 0.0;
  double rho = 1.0;
  for (int ic = 0; ic < nc; ic++) {
    /* Compute the modified local velocity */
    normUT = sqrt(w[ic * 3 + 1] * w[ic * 3 + 1] + w[ic * 3 + 2] * w[ic * 3 + 2]);
    worktmp[0] = w[ic * 3] + mu[ic] * normUT - rho * (z[ic * 3]);
    worktmp[1] = w[ic * 3 + 1] - rho * z[ic * 3 + 1];
    worktmp[2] = w[ic * 3 + 2] - rho * z[ic * 3 + 2];
    invmu = 1.0 / mu[ic];
    soc_projection_3d(worktmp, invmu);
    normUT = sqrt(worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
    worktmp[0] = w[ic * 3] - (worktmp[0] - mu[ic] * normUT);
    worktmp[1] = w[ic * 3 + 1] - worktmp[1];
    worktmp[2] = w[ic * 3 + 2] - worktmp[2];
    *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }
  *error = sqrt(*error);

  /* Computes error */
  double norm_q = cblas_dnrm2(n, problem->q, incx);
  *error = *error / (norm_q + 1.0);
  if (*error > tolerance) {
    /*      if (verbose > 0) printf(" Numerics - fc3d_compute_error_velocity failed: error = %g
     * > tolerance = %g.\n",*error, tolerance); */
    return 1;
  } else
    return 0;
}

void fc3d_Tresca_unitary_compute_and_add_error(double z[3], double w[3], double R,
                                               double* error, double* worktmp) {
  /* natural map z - P_C(z - w), C = R+ x B(0, R) the Tresca cylinder */
  (void)worktmp;
  ConvexSet K = convex_set(CONVEX_SET_POSITIVE_CYLINDER, 3, 1, &R, NULL);
  *error += convex_set_natural_map_block(&K, 0, z, w, NULL);
}
int fc3d_Tresca_compute_error(FrictionContactProblem* problem, double* z, double* w,
                              double tolerance, SolverOptions* options, double* error) {
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    return numerics_error("fc3d_Tresca_compute_error",
                          "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  double R;
  cblas_dcopy(n, problem->q, incx, w, incy);  // w <-q
  // Compute the current velocity
  /* NM_gemv(1.0, problem->M, z, 1.0, w); */
  NM_prod_mv_3x3(n, n, problem->M, z, w);
  *error = 0.;
  int ic, ic3;
  double worktmp[3];
  for (ic = 0, ic3 = 0; ic < nc; ic++, ic3 += 3) {
    R = (options->dWork[ic]);
    fc3d_Tresca_unitary_compute_and_add_error(z + ic3, w + ic3, R, error, worktmp);
  }
  *error = sqrt(*error);

  /* Computes error */
  double norm = cblas_dnrm2(n, problem->q, incx);
  DEBUG_PRINTF("norm = %12.8e\n", norm);
  if (fabs(norm) > DBL_EPSILON) *error /= norm;

  if (*error > tolerance) {
    /* if (verbose > 0) printf(" Numerics - fc3d_Tresca_compute_error failed: error = %g >
     * tolerance = %g.\n",*error, tolerance);  */
    return 1;
  } else
    return 0;
}
