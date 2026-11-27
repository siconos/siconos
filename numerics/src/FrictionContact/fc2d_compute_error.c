// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "fc2d_compute_error.h"

#include <assert.h>
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, sqrt
#include <stdio.h>   // for printf
#include <stdlib.h>  // for malloc, free

#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "FrictionContact_options.h"
#include "NumericsMatrix.h"     // for NM_gemv
#include "SiconosBlas.h"        // for cblas_dcopy, cblas_dnrm2, cblas_dscal
#include "SolverOptions.h"      // for SolverOptions
#include "convex_set.h"         // for ConvexSet
#include "error_computation.h"  // for error_compute_alloc_pair
#include "natural_map_error.h"  // for natural_map_error
#include "numerics_errors.h"
#include "numerics_verbose.h"  // for numerics_error, verbose

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

void fc2d_unitary_compute_and_add_error(const double* restrict r, const double* restrict u,
                                        double mu, double* restrict error, double* worktmp) {
  /* natural map r - P_K(r - u_mod), u_mod = u + mu |u_T| e_N the modified velocity */
  worktmp[0] = u[0] + mu * fabs(u[1]);
  worktmp[1] = u[1];
  ConvexSet K = convex_set(CONVEX_SET_SOC, 2, 1, &mu, NULL);
  *error += convex_set_natural_map_block(&K, 0, r, worktmp, NULL);
}

void fc2d_compute_modified_local_velocity(int nc, const double* u, const double* mu,
                                          double* u_modified, double* work) {
  int alloc_work = 0;
  if (work == NULL) {
    work = (double*)malloc(nc * sizeof(double));
    alloc_work = 1;
  }

  for (int i = 0; i < nc; ++i) {
    work[i] = fabs(u[2 * i + 1]);
  }

  for (int i = 0; i < nc; ++i) {
    const int pos = 2 * i;
    u_modified[pos] = u[pos] + mu[i] * work[i];
    u_modified[pos + 1] = u[pos + 1];
  }

  if (alloc_work) free(work);
}

double fc2d_compute_complementarity_error(int nc, const double* r, const double* u,
                                          const double* mu, int norm, int scope, int reduction,
                                          double alpha, double* work, double* tmp) {
  /* modified local velocity u_mod in tmp, then the scaled natural-map error on the stack of
     2D friction cones (natural_map_error.h) */
  fc2d_compute_modified_local_velocity(nc, u, mu, tmp, work);
  NaturalMapErrorScaling s = {norm, scope, reduction, alpha};
  ConvexSet K = convex_set(CONVEX_SET_SOC, 2, nc, mu, NULL);
  return natural_map_error(&s, &K, r, tmp);
}

int fc2d_compute_error(FrictionContactProblem* problem, double* z, double* w, double tolerance,
                       SolverOptions* options, double* error) {
  DEBUG_BEGIN("fc2d_compute_error(...)\n");
  /* Input validation using standardized macros */
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(w);
  CHECK_NULL(error);
  CHECK_MATRIX(problem->M);
  CHECK_NULL(problem->q);
  CHECK_NULL(problem->mu);

  int nc = problem->numberOfContacts;
  int n = nc * 2;

  double* mu = problem->mu;

  cblas_dcopy(n, problem->q, 1, w, 1);  // w <-q
  NM_gemv(1.0, problem->M, z, 1.0,
          w);  // w = W z +q is assumed to be compute at machine accuracy.

  if (nc == 0) {
    *error = 0.0;
    return 0;
  }

  const NaturalMapErrorScaling defaults = {SICONOS_FRICTION_2D_SCALING_EUCLIDEAN,
                                           SICONOS_FRICTION_2D_SCALING_SCOPE_PER_CONE,
                                           SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX, 1.0};
  NaturalMapErrorScaling sc = natural_map_error_scaling(options, defaults);
  int norm = sc.norm;

  double* comp_work = NULL;
  double* comp_tmp = NULL;
  if (error_compute_alloc_pair(nc, n, &comp_work, &comp_tmp)) {
    return numerics_error("fc2d_compute_error", "bad alloc");
  }
  *error = fc2d_compute_complementarity_error(nc, z, w, mu, sc.norm, sc.scope, sc.reduction,
                                              sc.alpha, comp_work, comp_tmp);
  free(comp_work);
  free(comp_tmp);

  /* Explicit no-scaling request: preserve the historical FC2D error
     definition error / max(norm_r, norm_w, norm_q). */
  if (norm == SICONOS_FRICTION_2D_SCALING_NO) {
    double norm_r = cblas_dnrm2(n, z, 1);
    double norm_w = cblas_dnrm2(n, w, 1);
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double relative_scaling = fmax(norm_r, fmax(norm_w, norm_q));
    if (fabs(relative_scaling) > DBL_EPSILON) *error /= relative_scaling;
  }

  DEBUG_PRINTF("error=%e\n", *error);
  DEBUG_END("fc2d_compute_error(...)\n");
  if (*error > tolerance) {
    numerics_printf_verbose(1, "-- FC2D - warning. fc2d_compute_error = %g > %g (= tolerance)",
                            *error, tolerance);
    return 1;
  } else
    return 0;
}
