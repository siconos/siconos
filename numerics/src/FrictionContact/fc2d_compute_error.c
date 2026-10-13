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
#include "fc2d_compute_error.h"

#include <assert.h>
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, sqrt
#include <stdio.h>   // for printf
#include <stdlib.h>  // for malloc, free

#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "FrictionContact_options.h"
#include "NumericsMatrix.h"  // for NM_gemv
#include "SOC_scaling.h"     // for SOC_scaling_*
#include "SiconosBlas.h"     // for cblas_dcopy, cblas_dnrm2, cblas_dscal
#include "SolverOptions.h"   // for SolverOptions
#include "numerics_errors.h"
#include "numerics_verbose.h"  // for numerics_error, verbose

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

void fc2d_unitary_compute_and_add_error(const double* restrict r, const double* restrict u,
                                        double mu, double* restrict error, double* worktmp) {
  /* Compute the modified local velocity */
  worktmp[0] = r[0] - u[0] - mu * fabs(u[1]);
  worktmp[1] = r[1] - u[1];
  /* projection */
  double normT = fabs(worktmp[1]);
  if (mu * normT <= -worktmp[0]) {
    worktmp[0] = 0.;
    worktmp[1] = 0.;
  } else if (normT > mu * worktmp[0]) {
    worktmp[0] = (mu * normT + worktmp[0]) / (mu * mu + 1.0);
    worktmp[1] = mu * worktmp[0] * SGN(worktmp[1]);
  }

  worktmp[0] = r[0] - worktmp[0];
  worktmp[1] = r[1] - worktmp[1];
  *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1];
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

static double safeguard_scale(double s) { return (s > 1.0) ? s : 1.0; }

/* Single-contact (2D cone) norm helpers. */
static inline double fc2d_unitary_norm(int norm, const double x[2], double mu, double alpha) {
  const double t = fabs(x[1]);
  switch (norm) {
    case SICONOS_FRICTION_2D_SCALING_EUCLIDEAN:
      return sqrt(x[0] * x[0] + t * t);
    case SICONOS_FRICTION_2D_SCALING_CONE:
      return fabs(x[0]) + t;
    case SICONOS_FRICTION_2D_SCALING_CONE_MU:
      return fabs(x[0]) + mu * t;
    case SICONOS_FRICTION_2D_SCALING_CONE_WEIGHTED:
      if (alpha < 0.0) alpha = 0.0;
      return sqrt(x[0] * x[0] + alpha * t * t);
    default:
      return 1.0;
  }
}

static inline double fc2d_unitary_scale(int norm, int reduction, const double r[2],
                                        const double u_mod[2], double mu, double alpha) {
  const double nr = fc2d_unitary_norm(norm, r, mu, alpha);
  const double nu = fc2d_unitary_norm(norm, u_mod, mu, alpha);

  switch (reduction) {
    case SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX: {
      const double sr = safeguard_scale(nr);
      const double su = safeguard_scale(nu);
      return (sr < su) ? sr : su;
    }
    case SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN:
      return sqrt(safeguard_scale(nr) * safeguard_scale(nu));
    case SICONOS_FRICTION_2D_SCALING_REDUCTION_MIXED:
      return nr + nu;
    default:
      return 1.0;
  }
}

/* Compute the scale for a stack of m 2D cones, using the selected norm and
   reduction.  Safeguarding (max(1, σ)) is applied to the final scale. */
static double fc2d_compute_vector_scale(int norm, int reduction, const double* x,
                                        const double* y, int m, const double* mu,
                                        double alpha) {
  double s = 1.0;

  switch (norm) {
    case SICONOS_FRICTION_2D_SCALING_EUCLIDEAN:
      switch (reduction) {
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_euclidean(x, y, m, 2);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_euclidean(x, y, m, 2);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_euclidean(x, y, m, 2);
          break;
      }
      break;

    case SICONOS_FRICTION_2D_SCALING_CONE:
      switch (reduction) {
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone(x, y, m, 2);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone(x, y, m, 2);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone(x, y, m, 2);
          break;
      }
      break;

    case SICONOS_FRICTION_2D_SCALING_CONE_MU:
      switch (reduction) {
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone_mu(x, y, m, mu, 2);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone_mu(x, y, m, mu, 2);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone_mu(x, y, m, mu, 2);
          break;
      }
      break;

    case SICONOS_FRICTION_2D_SCALING_CONE_WEIGHTED:
      switch (reduction) {
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone_weighted(x, y, m, 2, alpha);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone_weighted(x, y, m, 2, alpha);
          break;
        case SICONOS_FRICTION_2D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone_weighted(x, y, m, 2, alpha);
          break;
      }
      break;

    default:
      s = 1.0;
      break;
  }

  return safeguard_scale(s);
}

/* Project a 2D vector w onto the cone K = { (wN,wT) | |wT| <= mu*wN }.
   The result is written into p. */
static inline void fc2d_project_on_cone(const double w[2], double mu, double p[2]) {
  const double normT = fabs(w[1]);
  if (mu * normT <= -w[0]) {
    p[0] = 0.0;
    p[1] = 0.0;
  } else if (normT <= mu * w[0]) {
    p[0] = w[0];
    p[1] = w[1];
  } else {
    const double mu2 = mu * mu + 1.0;
    p[0] = (mu * normT + w[0]) / mu2;
    p[1] = mu * p[0] * SGN(w[1]);
  }
}

double fc2d_compute_complementarity_error(int nc, const double* r, const double* u,
                                          const double* mu, int norm, int scope, int reduction,
                                          double alpha, double* work, double* tmp) {
  const int m = 2 * nc;
  const int do_scale = (norm != SICONOS_FRICTION_2D_SCALING_NO);

  /* Build the modified local velocity u_mod once in tmp. */
  fc2d_compute_modified_local_velocity(nc, u, mu, tmp, work);

  if (!do_scale || scope == SICONOS_FRICTION_2D_SCALING_SCOPE_PER_VECTOR) {
    double global_scale = 1.0;
    if (do_scale) {
      global_scale = fc2d_compute_vector_scale(norm, reduction, r, tmp, nc, mu, alpha);
    }

    /* tmp <- r - u_mod, then project each contact. */
    for (int i = 0; i < m; ++i) {
      tmp[i] = r[i] - tmp[i];
    }
    for (int i = 0; i < nc; ++i) {
      const int pos = 2 * i;
      double p[2];
      fc2d_project_on_cone(&tmp[pos], mu[i], p);
      tmp[pos] = p[0];
      tmp[pos + 1] = p[1];
    }

    /* tmp <- unscaled natural-map residual, then scale the whole vector at once. */
    for (int i = 0; i < m; ++i) {
      tmp[i] = r[i] - tmp[i];
    }
    cblas_dscal(m, 1.0 / global_scale, tmp, 1);
  } else {
    /* Per-contact scaling: compute scale, projection and residual in one pass.
       The modified local velocity is read from tmp before it is overwritten. */
    for (int i = 0; i < nc; ++i) {
      const int pos = 2 * i;
      const double* r_i = &r[pos];
      const double* u_mod_i = &tmp[pos];

      const double s =
          safeguard_scale(fc2d_unitary_scale(norm, reduction, r_i, u_mod_i, mu[i], alpha));
      const double inv = 1.0 / s;

      const double w[2] = {r_i[0] - u_mod_i[0], r_i[1] - u_mod_i[1]};
      double p[2];
      fc2d_project_on_cone(w, mu[i], p);

      tmp[pos] = (r_i[0] - p[0]) * inv;
      tmp[pos + 1] = (r_i[1] - p[1]) * inv;
    }
  }

  return cblas_dnrm2(m, tmp, 1);
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

  int norm = SICONOS_FRICTION_2D_SCALING_EUCLIDEAN;
  int scope = SICONOS_FRICTION_2D_SCALING_SCOPE_PER_CONE;
  int reduction = SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX;
  double alpha = 1.0;
  if (options != NULL) {
    norm = options->iparam[SICONOS_FRICTION_2D_SCALING_NORM];
    scope = options->iparam[SICONOS_FRICTION_2D_SCALING_SCOPE];
    reduction = options->iparam[SICONOS_FRICTION_2D_SCALING_REDUCTION];
    alpha = options->dparam[SICONOS_FRICTION_2D_SCALING_ALPHA];
  }

  double* comp_work = (double*)malloc(nc * sizeof(double));
  double* comp_tmp = (double*)malloc(n * sizeof(double));
  if (comp_work == NULL || comp_tmp == NULL) {
    free(comp_work);
    free(comp_tmp);
    return numerics_error("fc2d_compute_error", "bad alloc");
  }
  *error = fc2d_compute_complementarity_error(nc, z, w, mu, norm, scope, reduction, alpha,
                                              comp_work, comp_tmp);
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
