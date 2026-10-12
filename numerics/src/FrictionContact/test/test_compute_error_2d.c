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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "SiconosBlas.h"             // for cblas_dnrm2
#include "SOC_scaling.h"             // for SOC_scaling_*
#include "fc2d_compute_error.h"
#include "FrictionContact_options.h"

#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

static double safeguard_scale(double s) { return (s > 1.0) ? s : 1.0; }

static double scalar_norm(int norm, const double x[2], double mu, double alpha) {
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

static double scalar_scale(int norm, int reduction, const double r[2],
                           const double u_mod[2], double mu, double alpha) {
  const double nr = scalar_norm(norm, r, mu, alpha);
  const double nu = scalar_norm(norm, u_mod, mu, alpha);
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

static void project_2d(const double w[2], double mu, double p[2]) {
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

static double scalar_complementarity_error(int nc, const double* r, const double* u,
                                           const double* mu, int norm, int scope,
                                           int reduction, double alpha) {
  double err_sq = 0.0;
  double global_scale = 1.0;

  if (norm != SICONOS_FRICTION_2D_SCALING_NO &&
      scope == SICONOS_FRICTION_2D_SCALING_SCOPE_PER_VECTOR) {
    double* u_mod = (double*)malloc(2 * nc * sizeof(double));
    double* work = (double*)malloc(nc * sizeof(double));
    fc2d_compute_modified_local_velocity(nc, u, mu, u_mod, work);
    switch (norm) {
      case SICONOS_FRICTION_2D_SCALING_EUCLIDEAN:
        global_scale = (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX)
                           ? SOC_scaling_euclidean(r, u_mod, nc, 2)
                           : (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN)
                                 ? SOC_geom_scaling_euclidean(r, u_mod, nc, 2)
                                 : SOC_mixed_scaling_euclidean(r, u_mod, nc, 2);
        break;
      case SICONOS_FRICTION_2D_SCALING_CONE:
        global_scale = (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX)
                           ? SOC_scaling_cone(r, u_mod, nc, 2)
                           : (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN)
                                 ? SOC_geom_scaling_cone(r, u_mod, nc, 2)
                                 : SOC_mixed_scaling_cone(r, u_mod, nc, 2);
        break;
      case SICONOS_FRICTION_2D_SCALING_CONE_MU:
        global_scale = (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX)
                           ? SOC_scaling_cone_mu(r, u_mod, nc, mu, 2)
                           : (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN)
                                 ? SOC_geom_scaling_cone_mu(r, u_mod, nc, mu, 2)
                                 : SOC_mixed_scaling_cone_mu(r, u_mod, nc, mu, 2);
        break;
      case SICONOS_FRICTION_2D_SCALING_CONE_WEIGHTED:
        global_scale = (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX)
                           ? SOC_scaling_cone_weighted(r, u_mod, nc, 2, alpha)
                           : (reduction == SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN)
                                 ? SOC_geom_scaling_cone_weighted(r, u_mod, nc, 2, alpha)
                                 : SOC_mixed_scaling_cone_weighted(r, u_mod, nc, 2, alpha);
        break;
    }
    global_scale = safeguard_scale(global_scale);
    free(u_mod);
    free(work);
  }

  for (int i = 0; i < nc; ++i) {
    const int pos = 2 * i;
    const double u_t = fabs(u[pos + 1]);
    const double u_mod[2] = {u[pos] + mu[i] * u_t, u[pos + 1]};

    double s = global_scale;
    if (norm != SICONOS_FRICTION_2D_SCALING_NO &&
        scope == SICONOS_FRICTION_2D_SCALING_SCOPE_PER_CONE) {
      s = safeguard_scale(scalar_scale(norm, reduction, &r[pos], u_mod, mu[i], alpha));
    }

    const double w[2] = {r[pos] - u_mod[0], r[pos + 1] - u_mod[1]};
    double p[2];
    project_2d(w, mu[i], p);

    const double inv = 1.0 / s;
    const double res0 = (r[pos] - p[0]) * inv;
    const double res1 = (r[pos + 1] - p[1]) * inv;
    err_sq += res0 * res0 + res1 * res1;
  }

  return sqrt(err_sq);
}

static double rel_diff(double a, double b) {
  if (a == 0.0 && b == 0.0) return 0.0;
  return fabs(a - b) / (fabs(a) + fabs(b) + 1.0);
}

int main(void) {
  int info = 0;
  const int nc = 8;
  const int n = 2 * nc;

  srand48(1234);
  double* r = (double*)malloc(n * sizeof(double));
  double* u = (double*)malloc(n * sizeof(double));
  double* mu = (double*)malloc(nc * sizeof(double));
  for (int i = 0; i < n; ++i) {
    /* Keep component magnitudes away from zero so that per-contact scales are
       not clipped by the safeguard(·) and scaling invariance holds exactly. */
    r[i] = (1.0 + drand48()) * (drand48() > 0.5 ? 1.0 : -1.0);
    u[i] = (1.0 + drand48()) * (drand48() > 0.5 ? 1.0 : -1.0);
  }
  for (int i = 0; i < nc; ++i) {
    mu[i] = 0.1 + 0.9 * drand48();
  }

  double* work = (double*)malloc(nc * sizeof(double));
  double* tmp = (double*)malloc(n * sizeof(double));

  int norms[] = {SICONOS_FRICTION_2D_SCALING_NO,       SICONOS_FRICTION_2D_SCALING_EUCLIDEAN,
                 SICONOS_FRICTION_2D_SCALING_CONE,    SICONOS_FRICTION_2D_SCALING_CONE_MU,
                 SICONOS_FRICTION_2D_SCALING_CONE_WEIGHTED};
  const char* norm_names[] = {"NO", "EUCLIDEAN", "CONE", "CONE_MU", "WEIGHTED"};
  const int nnorms = sizeof(norms) / sizeof(norms[0]);

  int scopes[] = {SICONOS_FRICTION_2D_SCALING_SCOPE_PER_CONE,
                  SICONOS_FRICTION_2D_SCALING_SCOPE_PER_VECTOR};
  const char* scope_names[] = {"PER_CONE", "PER_VECTOR"};
  const int nscopes = sizeof(scopes) / sizeof(scopes[0]);

  int reductions[] = {SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX,
                      SICONOS_FRICTION_2D_SCALING_REDUCTION_MEAN,
                      SICONOS_FRICTION_2D_SCALING_REDUCTION_MIXED};
  const char* reduction_names[] = {"MAX", "MEAN", "MIXED"};
  const int nreductions = sizeof(reductions) / sizeof(reductions[0]);

  for (int in = 0; in < nnorms; ++in) {
    for (int is = 0; is < nscopes; ++is) {
      for (int ir = 0; ir < nreductions; ++ir) {
        char label[64];
        snprintf(label, sizeof(label), "%s_%s_%s", norm_names[in], scope_names[is],
                 reduction_names[ir]);
        fprintf(stderr, "  scaling %s\n", label);

        double err_vec =
            fc2d_compute_complementarity_error(nc, r, u, mu, norms[in], scopes[is],
                                               reductions[ir], /*alpha=*/2.0, work, tmp);
        double err_ref = scalar_complementarity_error(nc, r, u, mu, norms[in], scopes[is],
                                                      reductions[ir], /*alpha=*/2.0);
        double rd = rel_diff(err_vec, err_ref);
        printf("fc2d %-40s rel diff = %e\n", label, rd);
        if (rd > 1e-12) {
          fprintf(stderr, "Mismatch for scaling %s!\n", label);
          info += 1;
        }

        /* Scaling invariance for norm != NO. */
        if (norms[in] != SICONOS_FRICTION_2D_SCALING_NO) {
          const double lambda = 3.7;
          double* r_scaled = (double*)malloc(n * sizeof(double));
          double* u_scaled = (double*)malloc(n * sizeof(double));
          for (int k = 0; k < n; ++k) {
            r_scaled[k] = lambda * r[k];
            u_scaled[k] = lambda * u[k];
          }
          double err_scaled = fc2d_compute_complementarity_error(nc, r_scaled, u_scaled, mu,
                                                                 norms[in], scopes[is],
                                                                 reductions[ir], 2.0, work, tmp);
          double rd_scale = rel_diff(err_vec, err_scaled);
          printf("fc2d %-40s scale invariance = %e\n", label, rd_scale);
          if (rd_scale > 1e-12) {
            fprintf(stderr, "Scaling invariance violated for %s!\n", label);
            info += 1;
          }
          free(r_scaled);
          free(u_scaled);
        }
      }
    }
  }

  free(r);
  free(u);
  free(mu);
  free(work);
  free(tmp);
  return info;
}
