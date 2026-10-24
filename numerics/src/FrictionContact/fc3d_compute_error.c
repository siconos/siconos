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
#include "SOC_scaling.h"  // for SOC_scaling_*
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "projectionOnCone.h"      // for projectionOnCone, projectionOnConeVector
#include "projectionOnCylinder.h"  // for projectionOnCylinder
#include "siconos_debug.h"         // for DEBUG_PRINTF, DEBUG_EXPR, DEBUG_...
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif
#include "SiconosBlas.h"  // for cblas_dcopy, cblas_dnrm2

void fc3d_unitary_compute_and_add_error(const double r[3], const double u[3], double mu,
                                        double* restrict error, double worktmp[3]) {
  // double normUT;
  // double worktmp[3];
  /* Compute the modified local velocity */
  /* worktmp[0] = r[0] - u[0] - mu *  hypot(u[1], u[2]); */
  worktmp[0] = r[0] - u[0] - mu * sqrt(u[1] * u[1] + u[2] * u[2]);
  worktmp[1] = r[1] - u[1];
  worktmp[2] = r[2] - u[2];
  projectionOnCone(worktmp, mu);
  worktmp[0] = r[0] - worktmp[0];
  worktmp[1] = r[1] - worktmp[1];
  worktmp[2] = r[2] - worktmp[2];
  *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
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

static double safeguard_scale(double s) { return (s > 1.0) ? s : 1.0; }

/* Helpers for the single-contact (3D cone) case.  These inline functions
   reproduce the per-variable scalings implemented in SOC_scaling.c for m=1
   without the overhead of a function call or BLAS on a 3-vector. */
static inline double fc3d_unitary_norm(int norm, const double x[3], double mu, double alpha) {
  const double t = sqrt(x[1] * x[1] + x[2] * x[2]);
  switch (norm) {
    case SICONOS_FRICTION_3D_SCALING_EUCLIDEAN:
      return sqrt(x[0] * x[0] + t * t);
    case SICONOS_FRICTION_3D_SCALING_CONE:
      return fabs(x[0]) + t;
    case SICONOS_FRICTION_3D_SCALING_CONE_MU:
      return fabs(x[0]) + mu * t;
    case SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED:
      if (alpha < 0.0) alpha = 0.0;
      return sqrt(x[0] * x[0] + alpha * t * t);
    default:
      return 1.0;
  }
}

static inline double fc3d_unitary_scale(int norm, int reduction, const double r[3],
                                        const double u_mod[3], double mu, double alpha) {
  const double nr = fc3d_unitary_norm(norm, r, mu, alpha);
  const double nu = fc3d_unitary_norm(norm, u_mod, mu, alpha);

  switch (reduction) {
    case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX: {
      const double sr = safeguard_scale(nr);
      const double su = safeguard_scale(nu);
      return (sr < su) ? sr : su;
    }
    case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
      return sqrt(safeguard_scale(nr) * safeguard_scale(nu));
    case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
      return nr + nu;
    default:
      return 1.0;
  }
}

/* Compute the scale for a stack of m cones of dimension dim, using the
   selected norm, scope and reduction.  Safeguarding (max(1, σ)) is applied
   to the final scale. */
static double fc3d_compute_vector_scale(int norm, int reduction, const double* x,
                                        const double* y, int m, int dim, const double* mu,
                                        double alpha) {
  double s = 1.0;

  switch (norm) {
    case SICONOS_FRICTION_3D_SCALING_EUCLIDEAN:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_euclidean(x, y, m, dim);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_euclidean(x, y, m, dim);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_euclidean(x, y, m, dim);
          break;
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone(x, y, m, dim);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone(x, y, m, dim);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone(x, y, m, dim);
          break;
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE_MU:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone_mu(x, y, m, mu, dim);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone_mu(x, y, m, mu, dim);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone_mu(x, y, m, mu, dim);
          break;
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone_weighted(x, y, m, dim, alpha);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone_weighted(x, y, m, dim, alpha);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone_weighted(x, y, m, dim, alpha);
          break;
      }
      break;

    default:
      s = 1.0;
      break;
  }

  return safeguard_scale(s);
}

double fc3d_compute_complementarity_error(int nc, const double* r, const double* u,
                                          const double* mu, int norm, int scope, int reduction,
                                          double alpha, double* work, double* tmp) {
  const int m = 3 * nc;
  const int do_scale = (norm != SICONOS_FRICTION_3D_SCALING_NO);

  /* Build the modified local velocity u_mod once in tmp.  It is reused by the
     whole-vector scaling path and by the per-cone scaling loop below. */
  fc3d_compute_modified_local_velocity(nc, u, mu, tmp, work);

  if (!do_scale || scope == SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR) {
    /* Whole-vector scaling: a single global scale is enough. */
    double global_scale = 1.0;
    if (do_scale) {
      global_scale = fc3d_compute_vector_scale(norm, reduction, r, tmp, nc, 3, mu, alpha);
    }

    /* tmp <- r - u_mod, then project to evaluate the unscaled natural map. */
    for (int i = 0; i < m; ++i) {
      tmp[i] = r[i] - tmp[i];
    }
    projectionOnConeVector(nc, tmp, mu, work);

    /* tmp <- unscaled natural-map residual, then scale the whole vector at once. */
    for (int i = 0; i < m; ++i) {
      tmp[i] = r[i] - tmp[i];
    }
    cblas_dscal(m, 1.0 / global_scale, tmp, 1);
  } else {
    /* Per-contact scaling: compute scale, projection and residual in one pass.
       The modified local velocity is read from tmp before it is overwritten. */
    for (int i = 0; i < nc; ++i) {
      const int pos = 3 * i;
      const double* r_i = &r[pos];
      const double* u_mod_i = &tmp[pos];

      const double s =
          safeguard_scale(fc3d_unitary_scale(norm, reduction, r_i, u_mod_i, mu[i], alpha));
      const double inv = 1.0 / s;

      /* w = r_i - u_mod_i, projection of w on the cone, residual = r_i - proj. */
      const double w0 = r_i[0] - u_mod_i[0];
      const double w1 = r_i[1] - u_mod_i[1];
      const double w2 = r_i[2] - u_mod_i[2];

      const double normTw = sqrt(w1 * w1 + w2 * w2);
      double p0, p1, p2;
      if (mu[i] * normTw <= -w0) {
        p0 = 0.0;
        p1 = 0.0;
        p2 = 0.0;
      } else if (normTw <= mu[i] * w0) {
        p0 = w0;
        p1 = w1;
        p2 = w2;
      } else {
        const double mu2 = mu[i] * mu[i] + 1.0;
        p0 = (mu[i] * normTw + w0) / mu2;
        const double scale = mu[i] * p0 / normTw;
        p1 = scale * w1;
        p2 = scale * w2;
      }

      tmp[pos] = (r_i[0] - p0) * inv;
      tmp[pos + 1] = (r_i[1] - p1) * inv;
      tmp[pos + 2] = (r_i[2] - p2) * inv;
    }
  }

  return cblas_dnrm2(m, tmp, 1);
}

void fc3d_unitary_compute_dual_and_add_error(double* restrict r, double* restrict u, double mu,
                                             double* restrict error, double* worktmp) {
  worktmp[0] = u[0] + mu * sqrt(u[1] * u[1] + u[2] * u[2]) - r[0];
  worktmp[1] = u[1] - r[1];
  worktmp[2] = u[2] - r[2];
  projectionOnDualCone(worktmp, mu);
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

  int norm_id = SICONOS_FRICTION_3D_SCALING_EUCLIDEAN;
  int scope = SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE;
  int reduction = SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX;
  double alpha = 1.0;
  if (options != NULL) {
    norm_id = options->iparam[SICONOS_FRICTION_3D_SCALING_NORM];
    scope = options->iparam[SICONOS_FRICTION_3D_SCALING_SCOPE];
    reduction = options->iparam[SICONOS_FRICTION_3D_SCALING_REDUCTION];
    alpha = options->dparam[SICONOS_FRICTION_3D_SCALING_ALPHA];
  }
  double* comp_work = NULL;
  double* comp_tmp = NULL;
  if (error_compute_alloc_pair(nc, n, &comp_work, &comp_tmp)) {
    return numerics_error("fc3d_compute_error", "bad alloc");
  }
  *error = fc3d_compute_complementarity_error(nc, z, w, mu, norm_id, scope, reduction, alpha,
                                              comp_work, comp_tmp);
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
    projectionOnCone(worktmp, invmu);
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
  /* Compute the modified local velocity */
  worktmp[0] = z[0] - w[0];
  worktmp[1] = z[1] - w[1];
  worktmp[2] = z[2] - w[2];
  projectionOnCylinder(worktmp, R);
  worktmp[0] = z[0] - worktmp[0];
  worktmp[1] = z[1] - worktmp[1];
  worktmp[2] = z[2] - worktmp[2];
  *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
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
