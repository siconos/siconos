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

#include "SOC_scaling.h"

#include <math.h>    // for sqrt, fabs, fmax
#include <stddef.h>  // for NULL

#include "SiconosBlas.h"  // for cblas_dnrm2, cblas_dasum, cblas_ddot

/* -------------------------------------------------------------------------- */
/* Single-vector norms                                                        */
/* -------------------------------------------------------------------------- */

double SOC_norm_euclidean(const double* x, int m, int dim) {
  if (m <= 0 || dim <= 0 || x == NULL) return 0.0;
  return cblas_dnrm2(m * dim, x, 1);
}

double SOC_norm_cone(const double* x, int m, int dim) {
  if (m <= 0 || dim <= 0 || x == NULL) return 0.0;
  /* Normal part: sum of |x[j*dim]|, j = 0 .. m-1. */
  double normal = cblas_dasum(m, x, dim);
  if (dim == 1) return normal;
  /* Tangential part: sum of ||tangential_j||_2. */
  double tangential = 0.0;
  for (int j = 0; j < m; ++j) {
    tangential += cblas_dnrm2(dim - 1, x + j * dim + 1, 1);
  }
  return normal + tangential;
}

double SOC_norm_cone_mu(const double* x, int m, const double* mu, int dim) {
  if (m <= 0 || dim <= 0 || x == NULL || mu == NULL) return 0.0;
  double normal = cblas_dasum(m, x, dim);
  if (dim == 1) return normal;
  double tangential = 0.0;
  for (int j = 0; j < m; ++j) {
    tangential += mu[j] * cblas_dnrm2(dim - 1, x + j * dim + 1, 1);
  }
  return normal + tangential;
}

double SOC_norm_cone_weighted(const double* x, int m, int dim, double alpha) {
  if (m <= 0 || dim <= 0 || x == NULL) return 0.0;
  if (alpha < 0.0) alpha = 0.0;
  /* Normal part: sum of x[j*dim]^2. */
  double normal_sq = cblas_ddot(m, x, dim, x, dim);
  if (dim == 1) return sqrt(normal_sq);
  /* Tangential part: sum of squared tangential components. */
  double tangential_sq = 0.0;
  for (int j = 0; j < m; ++j) {
    tangential_sq += cblas_ddot(dim - 1, x + j * dim + 1, 1, x + j * dim + 1, 1);
  }
  return sqrt(normal_sq + alpha * tangential_sq);
}

/* -------------------------------------------------------------------------- */
/* Mixed scalings                                                             */
/* -------------------------------------------------------------------------- */

double SOC_mixed_scaling_euclidean(const double* x, const double* y, int m, int dim) {
  return SOC_norm_euclidean(x, m, dim) + SOC_norm_euclidean(y, m, dim);
}

double SOC_mixed_scaling_cone(const double* x, const double* y, int m, int dim) {
  return SOC_norm_cone(x, m, dim) + SOC_norm_cone(y, m, dim);
}

double SOC_mixed_scaling_cone_mu(const double* x, const double* y, int m, const double* mu,
                                 int dim) {
  return SOC_norm_cone_mu(x, m, mu, dim) + SOC_norm_cone_mu(y, m, mu, dim);
}

double SOC_mixed_scaling_cone_weighted(const double* x, const double* y, int m, int dim,
                                       double alpha) {
  return SOC_norm_cone_weighted(x, m, dim, alpha) + SOC_norm_cone_weighted(y, m, dim, alpha);
}

/* -------------------------------------------------------------------------- */
/* Per-variable scalings                                                      */
/* -------------------------------------------------------------------------- */

static double safeguard_scale(double s) { return (s > 1.0) ? s : 1.0; }

double SOC_scaling_euclidean(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard_scale(SOC_norm_euclidean(x, m, dim));
  const double sy = safeguard_scale(SOC_norm_euclidean(y, m, dim));
  return (sx < sy) ? sx : sy;
}

double SOC_scaling_cone(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard_scale(SOC_norm_cone(x, m, dim));
  const double sy = safeguard_scale(SOC_norm_cone(y, m, dim));
  return (sx < sy) ? sx : sy;
}

double SOC_scaling_cone_mu(const double* x, const double* y, int m, const double* mu,
                           int dim) {
  const double sx = safeguard_scale(SOC_norm_cone_mu(x, m, mu, dim));
  const double sy = safeguard_scale(SOC_norm_cone_mu(y, m, mu, dim));
  return (sx < sy) ? sx : sy;
}

double SOC_scaling_cone_weighted(const double* x, const double* y, int m, int dim,
                                 double alpha) {
  const double sx = safeguard_scale(SOC_norm_cone_weighted(x, m, dim, alpha));
  const double sy = safeguard_scale(SOC_norm_cone_weighted(y, m, dim, alpha));
  return (sx < sy) ? sx : sy;
}

/* -------------------------------------------------------------------------- */
/* Geometric per-variable scalings                                            */
/* -------------------------------------------------------------------------- */

double SOC_geom_scaling_euclidean(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard_scale(SOC_norm_euclidean(x, m, dim));
  const double sy = safeguard_scale(SOC_norm_euclidean(y, m, dim));
  return sqrt(sx * sy);
}

double SOC_geom_scaling_cone(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard_scale(SOC_norm_cone(x, m, dim));
  const double sy = safeguard_scale(SOC_norm_cone(y, m, dim));
  return sqrt(sx * sy);
}

double SOC_geom_scaling_cone_mu(const double* x, const double* y, int m, const double* mu,
                                int dim) {
  const double sx = safeguard_scale(SOC_norm_cone_mu(x, m, mu, dim));
  const double sy = safeguard_scale(SOC_norm_cone_mu(y, m, mu, dim));
  return sqrt(sx * sy);
}

double SOC_geom_scaling_cone_weighted(const double* x, const double* y, int m, int dim,
                                      double alpha) {
  const double sx = safeguard_scale(SOC_norm_cone_weighted(x, m, dim, alpha));
  const double sy = safeguard_scale(SOC_norm_cone_weighted(y, m, dim, alpha));
  return sqrt(sx * sy);
}
