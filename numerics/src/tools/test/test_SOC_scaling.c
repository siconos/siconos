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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "SOC_scaling.h"

static double max_abs(double a, double b) { return (fabs(a) > fabs(b)) ? a : b; }

static double rel_diff(double a, double b) {
  double denom = max_abs(a, b);
  return (denom > DBL_EPSILON) ? fabs(a - b) / denom : fabs(a - b);
}

static double scalar_tangential_norm(const double* z, int dim) {
  double s = 0.0;
  for (int i = 1; i < dim; ++i) s += z[i] * z[i];
  return sqrt(s);
}

static double scalar_norm_euclidean(const double* x, int m, int dim) {
  double s = 0.0;
  const int n = m * dim;
  for (int i = 0; i < n; ++i) s += x[i] * x[i];
  return sqrt(s);
}

static double scalar_norm_cone(const double* x, int m, int dim) {
  double s = 0.0;
  for (int j = 0; j < m; ++j) {
    const double* xj = x + j * dim;
    s += fabs(xj[0]) + scalar_tangential_norm(xj, dim);
  }
  return s;
}

static double scalar_norm_cone_mu(const double* x, int m, const double* mu, int dim) {
  double s = 0.0;
  for (int j = 0; j < m; ++j) {
    const double* xj = x + j * dim;
    s += fabs(xj[0]) + mu[j] * scalar_tangential_norm(xj, dim);
  }
  return s;
}

static double scalar_norm_cone_weighted(const double* x, int m, int dim, double alpha) {
  double s = 0.0;
  for (int j = 0; j < m; ++j) {
    const double* xj = x + j * dim;
    const double nt = scalar_tangential_norm(xj, dim);
    s += xj[0] * xj[0] + alpha * nt * nt;
  }
  return sqrt(s);
}

static double scalar_mixed_scaling_euclidean(const double* x, const double* y, int m,
                                             int dim) {
  return scalar_norm_euclidean(x, m, dim) + scalar_norm_euclidean(y, m, dim);
}

static double scalar_mixed_scaling_cone(const double* x, const double* y, int m, int dim) {
  return scalar_norm_cone(x, m, dim) + scalar_norm_cone(y, m, dim);
}

static double scalar_mixed_scaling_cone_mu(const double* x, const double* y, int m,
                                           const double* mu, int dim) {
  return scalar_norm_cone_mu(x, m, mu, dim) + scalar_norm_cone_mu(y, m, mu, dim);
}

static double scalar_mixed_scaling_cone_weighted(const double* x, const double* y, int m,
                                                 int dim, double alpha) {
  return scalar_norm_cone_weighted(x, m, dim, alpha) +
         scalar_norm_cone_weighted(y, m, dim, alpha);
}

static double safeguard(double s) { return (s > 1.0) ? s : 1.0; }

static double scalar_scaling_euclidean(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard(scalar_norm_euclidean(x, m, dim));
  const double sy = safeguard(scalar_norm_euclidean(y, m, dim));
  return (sx < sy) ? sx : sy;
}

static double scalar_scaling_cone(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard(scalar_norm_cone(x, m, dim));
  const double sy = safeguard(scalar_norm_cone(y, m, dim));
  return (sx < sy) ? sx : sy;
}

static double scalar_scaling_cone_mu(const double* x, const double* y, int m, const double* mu,
                                     int dim) {
  const double sx = safeguard(scalar_norm_cone_mu(x, m, mu, dim));
  const double sy = safeguard(scalar_norm_cone_mu(y, m, mu, dim));
  return (sx < sy) ? sx : sy;
}

static double scalar_scaling_cone_weighted(const double* x, const double* y, int m, int dim,
                                           double alpha) {
  const double sx = safeguard(scalar_norm_cone_weighted(x, m, dim, alpha));
  const double sy = safeguard(scalar_norm_cone_weighted(y, m, dim, alpha));
  return (sx < sy) ? sx : sy;
}

static double scalar_geom_scaling_euclidean(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard(scalar_norm_euclidean(x, m, dim));
  const double sy = safeguard(scalar_norm_euclidean(y, m, dim));
  return sqrt(sx * sy);
}

static double scalar_geom_scaling_cone(const double* x, const double* y, int m, int dim) {
  const double sx = safeguard(scalar_norm_cone(x, m, dim));
  const double sy = safeguard(scalar_norm_cone(y, m, dim));
  return sqrt(sx * sy);
}

static double scalar_geom_scaling_cone_mu(const double* x, const double* y, int m,
                                          const double* mu, int dim) {
  const double sx = safeguard(scalar_norm_cone_mu(x, m, mu, dim));
  const double sy = safeguard(scalar_norm_cone_mu(y, m, mu, dim));
  return sqrt(sx * sy);
}

static double scalar_geom_scaling_cone_weighted(const double* x, const double* y, int m,
                                                int dim, double alpha) {
  const double sx = safeguard(scalar_norm_cone_weighted(x, m, dim, alpha));
  const double sy = safeguard(scalar_norm_cone_weighted(y, m, dim, alpha));
  return sqrt(sx * sy);
}

static void fill_random(int n, double* v) {
  for (int i = 0; i < n; ++i) v[i] = 2.0 * drand48() - 1.0;
}

static int test_norms(void) {
  int info = 0;
  const int m = 2;
  const int dim = 3;
  const int n = m * dim;
  const int ntests = 20;

  srand48(42);
  for (int t = 0; t < ntests; ++t) {
    double x[6];
    fill_random(n, x);

    double ref, val;
    ref = scalar_norm_euclidean(x, m, dim);
    val = SOC_norm_euclidean(x, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_norm_euclidean mismatch\n");
      info += 1;
    }

    ref = scalar_norm_cone(x, m, dim);
    val = SOC_norm_cone(x, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_norm_cone mismatch\n");
      info += 1;
    }

    double mu[2] = {0.1 + 0.9 * drand48(), 0.1 + 0.9 * drand48()};
    ref = scalar_norm_cone_mu(x, m, mu, dim);
    val = SOC_norm_cone_mu(x, m, mu, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_norm_cone_mu mismatch\n");
      info += 1;
    }

    double alpha = 0.1 + 2.0 * drand48();
    ref = scalar_norm_cone_weighted(x, m, dim, alpha);
    val = SOC_norm_cone_weighted(x, m, dim, alpha);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_norm_cone_weighted mismatch\n");
      info += 1;
    }
  }

  printf("test_norms: %s\n", info == 0 ? "passed" : "FAILED");
  return info;
}

static int test_mixed_scalings(void) {
  int info = 0;
  const int m = 2;
  const int dim = 3;
  const int n = m * dim;
  const int ntests = 20;

  srand48(42);
  for (int t = 0; t < ntests; ++t) {
    double x[6], y[6];
    fill_random(n, x);
    fill_random(n, y);

    double ref, val;
    ref = scalar_mixed_scaling_euclidean(x, y, m, dim);
    val = SOC_mixed_scaling_euclidean(x, y, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_mixed_scaling_euclidean mismatch\n");
      info += 1;
    }

    ref = scalar_mixed_scaling_cone(x, y, m, dim);
    val = SOC_mixed_scaling_cone(x, y, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_mixed_scaling_cone mismatch\n");
      info += 1;
    }

    double mu[2] = {0.1 + 0.9 * drand48(), 0.1 + 0.9 * drand48()};
    ref = scalar_mixed_scaling_cone_mu(x, y, m, mu, dim);
    val = SOC_mixed_scaling_cone_mu(x, y, m, mu, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_mixed_scaling_cone_mu mismatch\n");
      info += 1;
    }

    double alpha = 0.1 + 2.0 * drand48();
    ref = scalar_mixed_scaling_cone_weighted(x, y, m, dim, alpha);
    val = SOC_mixed_scaling_cone_weighted(x, y, m, dim, alpha);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_mixed_scaling_cone_weighted mismatch\n");
      info += 1;
    }
  }

  printf("test_mixed_scalings: %s\n", info == 0 ? "passed" : "FAILED");
  return info;
}

static int test_per_variable_scalings(void) {
  int info = 0;
  const int m = 2;
  const int dim = 3;
  const int n = m * dim;
  const int ntests = 20;

  srand48(42);
  for (int t = 0; t < ntests; ++t) {
    double x[6], y[6];
    fill_random(n, x);
    fill_random(n, y);

    double ref, val;
    ref = scalar_scaling_euclidean(x, y, m, dim);
    val = SOC_scaling_euclidean(x, y, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_scaling_euclidean mismatch\n");
      info += 1;
    }

    ref = scalar_scaling_cone(x, y, m, dim);
    val = SOC_scaling_cone(x, y, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_scaling_cone mismatch\n");
      info += 1;
    }

    double mu[2] = {0.1 + 0.9 * drand48(), 0.1 + 0.9 * drand48()};
    ref = scalar_scaling_cone_mu(x, y, m, mu, dim);
    val = SOC_scaling_cone_mu(x, y, m, mu, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_scaling_cone_mu mismatch\n");
      info += 1;
    }

    double alpha = 0.1 + 2.0 * drand48();
    ref = scalar_scaling_cone_weighted(x, y, m, dim, alpha);
    val = SOC_scaling_cone_weighted(x, y, m, dim, alpha);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_scaling_cone_weighted mismatch\n");
      info += 1;
    }
  }

  printf("test_per_variable_scalings: %s\n", info == 0 ? "passed" : "FAILED");
  return info;
}

static int test_geom_scalings(void) {
  int info = 0;
  const int m = 2;
  const int dim = 3;
  const int n = m * dim;
  const int ntests = 20;

  srand48(42);
  for (int t = 0; t < ntests; ++t) {
    double x[6], y[6];
    fill_random(n, x);
    fill_random(n, y);

    double ref, val;
    ref = scalar_geom_scaling_euclidean(x, y, m, dim);
    val = SOC_geom_scaling_euclidean(x, y, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_geom_scaling_euclidean mismatch\n");
      info += 1;
    }

    ref = scalar_geom_scaling_cone(x, y, m, dim);
    val = SOC_geom_scaling_cone(x, y, m, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_geom_scaling_cone mismatch\n");
      info += 1;
    }

    double mu[2] = {0.1 + 0.9 * drand48(), 0.1 + 0.9 * drand48()};
    ref = scalar_geom_scaling_cone_mu(x, y, m, mu, dim);
    val = SOC_geom_scaling_cone_mu(x, y, m, mu, dim);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_geom_scaling_cone_mu mismatch\n");
      info += 1;
    }

    double alpha = 0.1 + 2.0 * drand48();
    ref = scalar_geom_scaling_cone_weighted(x, y, m, dim, alpha);
    val = SOC_geom_scaling_cone_weighted(x, y, m, dim, alpha);
    if (rel_diff(ref, val) > 1e-14) {
      fprintf(stderr, "SOC_geom_scaling_cone_weighted mismatch\n");
      info += 1;
    }
  }

  printf("test_geom_scalings: %s\n", info == 0 ? "passed" : "FAILED");
  return info;
}

static int check_homogeneity_one(const char* name, double s0, double sl, double lambda) {
  if (rel_diff(sl, lambda * s0) > 1e-14) {
    fprintf(stderr, "%s not homogeneous for lambda=%e\n", name, lambda);
    return 1;
  }
  return 0;
}

static int test_homogeneity(void) {
  int info = 0;
  const int m = 2;
  const int dim = 3;
  const int n = m * dim;
  double x[6], y[6];

  srand48(123);
  fill_random(n, x);
  fill_random(n, y);
  double mu[2] = {0.5, 0.3};
  double alpha = 2.0;

  double lambdas[] = {1e-8, 1e-4, 1.0, 1e4, 1e8};
  int nlam = sizeof(lambdas) / sizeof(lambdas[0]);

  for (int i = 0; i < nlam; ++i) {
    double lambda = lambdas[i];
    double xs[6], ys[6];
    for (int k = 0; k < n; ++k) {
      xs[k] = lambda * x[k];
      ys[k] = lambda * y[k];
    }

    info += check_homogeneity_one("SOC_norm_euclidean", SOC_norm_euclidean(x, m, dim),
                                  SOC_norm_euclidean(xs, m, dim), lambda);
    info += check_homogeneity_one("SOC_norm_cone", SOC_norm_cone(x, m, dim),
                                  SOC_norm_cone(xs, m, dim), lambda);
    info += check_homogeneity_one("SOC_norm_cone_mu", SOC_norm_cone_mu(x, m, mu, dim),
                                  SOC_norm_cone_mu(xs, m, mu, dim), lambda);
    info += check_homogeneity_one("SOC_norm_cone_weighted",
                                  SOC_norm_cone_weighted(x, m, dim, alpha),
                                  SOC_norm_cone_weighted(xs, m, dim, alpha), lambda);

    info += check_homogeneity_one("SOC_mixed_scaling_euclidean",
                                  SOC_mixed_scaling_euclidean(x, y, m, dim),
                                  SOC_mixed_scaling_euclidean(xs, ys, m, dim), lambda);
    info +=
        check_homogeneity_one("SOC_mixed_scaling_cone", SOC_mixed_scaling_cone(x, y, m, dim),
                              SOC_mixed_scaling_cone(xs, ys, m, dim), lambda);
    info += check_homogeneity_one("SOC_mixed_scaling_cone_mu",
                                  SOC_mixed_scaling_cone_mu(x, y, m, mu, dim),
                                  SOC_mixed_scaling_cone_mu(xs, ys, m, mu, dim), lambda);
    info +=
        check_homogeneity_one("SOC_mixed_scaling_cone_weighted",
                              SOC_mixed_scaling_cone_weighted(x, y, m, dim, alpha),
                              SOC_mixed_scaling_cone_weighted(xs, ys, m, dim, alpha), lambda);

    /* Per-variable scalings are safeguarded (max(1, ·)) and therefore are
       homogeneous only when the underlying norms are larger than 1. */
    if (lambda >= 1.0) {
      info +=
          check_homogeneity_one("SOC_scaling_euclidean", SOC_scaling_euclidean(x, y, m, dim),
                                SOC_scaling_euclidean(xs, ys, m, dim), lambda);
      info += check_homogeneity_one("SOC_scaling_cone", SOC_scaling_cone(x, y, m, dim),
                                    SOC_scaling_cone(xs, ys, m, dim), lambda);
      info +=
          check_homogeneity_one("SOC_scaling_cone_mu", SOC_scaling_cone_mu(x, y, m, mu, dim),
                                SOC_scaling_cone_mu(xs, ys, m, mu, dim), lambda);
      info += check_homogeneity_one("SOC_scaling_cone_weighted",
                                    SOC_scaling_cone_weighted(x, y, m, dim, alpha),
                                    SOC_scaling_cone_weighted(xs, ys, m, dim, alpha), lambda);

      info += check_homogeneity_one("SOC_geom_scaling_euclidean",
                                    SOC_geom_scaling_euclidean(x, y, m, dim),
                                    SOC_geom_scaling_euclidean(xs, ys, m, dim), lambda);
      info +=
          check_homogeneity_one("SOC_geom_scaling_cone", SOC_geom_scaling_cone(x, y, m, dim),
                                SOC_geom_scaling_cone(xs, ys, m, dim), lambda);
      info += check_homogeneity_one("SOC_geom_scaling_cone_mu",
                                    SOC_geom_scaling_cone_mu(x, y, m, mu, dim),
                                    SOC_geom_scaling_cone_mu(xs, ys, m, mu, dim), lambda);
      info +=
          check_homogeneity_one("SOC_geom_scaling_cone_weighted",
                                SOC_geom_scaling_cone_weighted(x, y, m, dim, alpha),
                                SOC_geom_scaling_cone_weighted(xs, ys, m, dim, alpha), lambda);
    }
  }

  printf("test_homogeneity: %s\n", info == 0 ? "passed" : "FAILED");
  return info;
}

int main(void) {
  int info = 0;
  info += test_norms();
  info += test_mixed_scalings();
  info += test_per_variable_scalings();
  info += test_geom_scalings();
  info += test_homogeneity();
  return info;
}
