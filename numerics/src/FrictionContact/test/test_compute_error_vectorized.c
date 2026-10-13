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
#include <string.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "GlobalFrictionContactProblem.h"
#include "NumericsMatrix.h"
#include "SOC_scaling.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "fc3d_compute_error.h"
#include "gfc3d_compute_error.h"
#include "projectionOnCone.h"

static double max_abs(double a, double b) { return a > b ? a : b; }

static double rel_diff(double a, double b) {
  double denom = max_abs(fabs(a), fabs(b));
  return (denom > DBL_EPSILON) ? fabs(a - b) / denom : fabs(a - b);
}

static NumericsMatrix* identity_dense(int n) {
  double* data = (double*)calloc((size_t)n * n, sizeof(double));
  for (int i = 0; i < n; ++i) data[i * n + i] = 1.0;
  return NM_create_from_data(NM_DENSE, n, n, data);
}

static void fill_random(int n, double* v) {
  for (int i = 0; i < n; ++i) v[i] = 2.0 * drand48() - 1.0;
}

/* Compute the scale for a stack of m cones (dim = 3), combining a norm,
   a scope and a reduction.  Safeguarding max(1, s) is applied. */
static double scalar_scale(int norm, int reduction, const double* x, const double* y, int m,
                           const double* mu, double alpha) {
  double s = 1.0;

  switch (norm) {
    case SICONOS_FRICTION_3D_SCALING_EUCLIDEAN:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_euclidean(x, y, m, 3);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_euclidean(x, y, m, 3);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_euclidean(x, y, m, 3);
          break;
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone(x, y, m, 3);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone(x, y, m, 3);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone(x, y, m, 3);
          break;
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE_MU:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone_mu(x, y, m, mu, 3);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone_mu(x, y, m, mu, 3);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone_mu(x, y, m, mu, 3);
          break;
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          s = SOC_scaling_cone_weighted(x, y, m, 3, alpha);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          s = SOC_geom_scaling_cone_weighted(x, y, m, 3, alpha);
          break;
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          s = SOC_mixed_scaling_cone_weighted(x, y, m, 3, alpha);
          break;
      }
      break;

    default:
      s = 1.0;
      break;
  }

  return (s > 1.0) ? s : 1.0;
}

/* Scalar reference: build u_mod, evaluate the unscaled natural map, then
   divide each contact residual by the scale dictated by norm/scope/reduction. */
static double scalar_complementarity_error(int nc, const double* r, const double* u,
                                           const double* mu, int norm, int scope,
                                           int reduction, double alpha, double norm_q) {
  (void)norm_q; /* unused: kept for API compatibility */
  const int n = 3 * nc;
  double* u_mod = (double*)malloc(n * sizeof(double));
  double* work = (double*)malloc(nc * sizeof(double));
  fc3d_compute_modified_local_velocity(nc, u, mu, u_mod, work);

  /* Whole-vector scaling uses a single global denominator. */
  double global_scale = 1.0;
  if (norm != SICONOS_FRICTION_3D_SCALING_NO &&
      scope == SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR) {
    global_scale = scalar_scale(norm, reduction, r, u_mod, nc, mu, alpha);
  }

  /* Accumulate the squared scaled residual. */
  double err_sq = 0.0;
  double worktmp[3];
  for (int i = 0; i < nc; ++i) {
    worktmp[0] = r[3 * i] - u_mod[3 * i];
    worktmp[1] = r[3 * i + 1] - u_mod[3 * i + 1];
    worktmp[2] = r[3 * i + 2] - u_mod[3 * i + 2];
    projectionOnCone(worktmp, mu[i]);
    worktmp[0] = r[3 * i] - worktmp[0];
    worktmp[1] = r[3 * i + 1] - worktmp[1];
    worktmp[2] = r[3 * i + 2] - worktmp[2];

    double s = global_scale;
    if (norm != SICONOS_FRICTION_3D_SCALING_NO &&
        scope == SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE) {
      s = scalar_scale(norm, reduction, &r[3 * i], &u_mod[3 * i], 1, &mu[i], alpha);
    }

    const double inv = 1.0 / s;
    err_sq += (worktmp[0] * inv) * (worktmp[0] * inv) +
              (worktmp[1] * inv) * (worktmp[1] * inv) +
              (worktmp[2] * inv) * (worktmp[2] * inv);
  }

  free(u_mod);
  free(work);
  return sqrt(err_sq);
}

static int test_fc3d_compute_error_scaling(void) {
  int info = 0;
  const int nc = 6;
  const int n = 3 * nc;

  srand48(123);
  double* r = (double*)malloc(n * sizeof(double));
  double* u = (double*)malloc(n * sizeof(double));
  double* mu = (double*)malloc(nc * sizeof(double));
  fill_random(n, r);
  fill_random(n, u);
  for (int i = 0; i < nc; ++i) mu[i] = 0.1 + 0.9 * drand48();

  double* work = (double*)malloc(nc * sizeof(double));
  double* tmp = (double*)malloc(n * sizeof(double));
  double norm_q = cblas_dnrm2(n, r, 1);

  int norms[] = {SICONOS_FRICTION_3D_SCALING_NO, SICONOS_FRICTION_3D_SCALING_EUCLIDEAN,
                 SICONOS_FRICTION_3D_SCALING_CONE, SICONOS_FRICTION_3D_SCALING_CONE_MU,
                 SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED};
  const char* norm_names[] = {"NO", "EUCLIDEAN", "CONE", "CONE_MU", "WEIGHTED"};
  const int nnorms = sizeof(norms) / sizeof(norms[0]);

  int scopes[] = {SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
                  SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR};
  const char* scope_names[] = {"PER_CONE", "PER_VECTOR"};
  const int nscopes = sizeof(scopes) / sizeof(scopes[0]);

  int reductions[] = {SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX,
                      SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN,
                      SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED};
  const char* reduction_names[] = {"MAX", "MEAN", "MIXED"};
  const int nreductions = sizeof(reductions) / sizeof(reductions[0]);

  for (int in = 0; in < nnorms; ++in) {
    for (int is = 0; is < nscopes; ++is) {
      for (int ir = 0; ir < nreductions; ++ir) {
        char label[64];
        snprintf(label, sizeof(label), "%s_%s_%s", norm_names[in], scope_names[is],
                 reduction_names[ir]);
        fprintf(stderr, "  scaling %s\n", label);
        double err_vec = fc3d_compute_complementarity_error(
            nc, r, u, mu, norms[in], scopes[is], reductions[ir], /*alpha=*/2.0, work, tmp);
        double err_ref =
            scalar_complementarity_error(nc, r, u, mu, norms[in], scopes[is], reductions[ir],
                                         /*alpha=*/2.0, norm_q);
        double rd = rel_diff(err_vec, err_ref);
        printf("fc3d %-40s rel diff = %e\n", label, rd);
        if (rd > 1e-12) {
          fprintf(stderr, "Mismatch for scaling %s!\n", label);
          info += 1;
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

static int test_fc3d_compute_error(void) {
  int info = 0;
  const int nc = 5;
  const int n = 3 * nc;

  FrictionContactProblem* problem = frictionContactProblem_new();
  problem->numberOfContacts = nc;
  problem->mu = (double*)malloc(nc * sizeof(double));
  problem->q = (double*)malloc(n * sizeof(double));
  problem->M = identity_dense(n);

  srand48(42);
  fill_random(n, problem->q);
  for (int i = 0; i < nc; ++i) problem->mu[i] = 0.1 + 0.9 * drand48();

  double* z = (double*)malloc(n * sizeof(double));
  double* w = (double*)malloc(n * sizeof(double));
  fill_random(n, z);

  /* w = M z + q */
  memcpy(w, problem->q, n * sizeof(double));
  NM_prod_mv_3x3(n, n, problem->M, z, w);

  SolverOptions options;
  memset(&options, 0, sizeof(options));
  options.iparam = (int*)calloc(OPTIONS_PARAM_SIZE, sizeof(int));
  options.dparam = (double*)calloc(OPTIONS_PARAM_SIZE, sizeof(double));
  options.iparam[SICONOS_FRICTION_3D_SCALING_NORM] = SICONOS_FRICTION_3D_SCALING_EUCLIDEAN;
  options.iparam[SICONOS_FRICTION_3D_SCALING_SCOPE] =
      SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE;
  options.iparam[SICONOS_FRICTION_3D_SCALING_REDUCTION] =
      SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX;
  options.dparam[SICONOS_FRICTION_3D_SCALING_ALPHA] = 1.0;

  double error = 0.0;
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  int err_info = fc3d_compute_error(problem, z, w, 1e-10, &options, &error);

  double error_ref = scalar_complementarity_error(
      nc, z, w, problem->mu, SICONOS_FRICTION_3D_SCALING_EUCLIDEAN,
      SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE, SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX,
      1.0, norm_q);
  double rd = rel_diff(error, error_ref);
  printf("fc3d_compute_error (EUCLIDEAN_PER_CONE_MAX) rel diff = %e (info=%d)\n", rd,
         err_info);
  if (rd > 1e-12) {
    fprintf(stderr, "fc3d_compute_error differs from scalar reference!\n");
    info += 1;
  }

  free(options.iparam);
  free(options.dparam);
  free(z);
  free(w);
  free(problem->mu);
  free(problem->q);
  free(problem->M->matrix0);
  problem->M->matrix0 = NULL;
  NM_free(problem->M);
  free(problem);

  return info;
}

static int test_gfc3d_compute_error(void) {
  int info = 0;
  const int nc = 4;
  const int m = 3 * nc;
  const int n = m;

  GlobalFrictionContactProblem* problem = globalFrictionContactProblem_new();
  problem->dimension = 3;
  problem->numberOfContacts = nc;
  problem->mu = (double*)malloc(nc * sizeof(double));
  problem->q = (double*)malloc(n * sizeof(double));
  problem->b = (double*)malloc(m * sizeof(double));
  problem->M = identity_dense(n);
  problem->H = identity_dense(m);

  srand48(77);
  fill_random(n, problem->q);
  fill_random(m, problem->b);
  for (int i = 0; i < nc; ++i) problem->mu[i] = 0.1 + 0.9 * drand48();

  double* reaction = (double*)malloc(m * sizeof(double));
  double* velocity = (double*)malloc(m * sizeof(double));
  double* globalVelocity = (double*)malloc(n * sizeof(double));
  fill_random(m, reaction);
  fill_random(n, globalVelocity);

  SolverOptions* options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
  if (!options) {
    fprintf(stderr, "Failed to create solver options\n");
    return 1;
  }
  options->iparam[SICONOS_FRICTION_3D_SCALING_NORM] = SICONOS_FRICTION_3D_SCALING_CONE_MU;
  options->iparam[SICONOS_FRICTION_3D_SCALING_SCOPE] =
      SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE;
  options->iparam[SICONOS_FRICTION_3D_SCALING_REDUCTION] =
      SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX;
  options->dparam[SICONOS_FRICTION_3D_SCALING_ALPHA] = 1.0;

  /* Reference dynamics error. */
  double* velocity_ref = (double*)malloc(m * sizeof(double));
  memcpy(velocity_ref, problem->b, m * sizeof(double));
  NM_tgemv(1, problem->H, globalVelocity, 1.0, velocity_ref);

  double* tmp = (double*)malloc(n * sizeof(double));
  double* tmp_1 = (double*)malloc(n * sizeof(double));
  memcpy(tmp_1, problem->q, n * sizeof(double));
  NM_gemv(1.0, problem->H, reaction, 0.0, tmp);
  double norm_Hr = cblas_dnrm2(n, tmp, 1);
  cblas_daxpy(n, 1.0, tmp, 1, tmp_1, 1);
  NM_gemv(-1.0, problem->M, globalVelocity, 0.0, tmp);
  double norm_Mv = cblas_dnrm2(n, tmp, 1);
  cblas_daxpy(n, 1.0, tmp, 1, tmp_1, 1);
  double err_dyn = cblas_dnrm2(n, tmp_1, 1);
  double relative_scaling = max_abs(cblas_dnrm2(n, problem->q, 1), max_abs(norm_Mv, norm_Hr));
  if (relative_scaling > DBL_EPSILON) err_dyn /= relative_scaling;

  double err_comp = scalar_complementarity_error(
      nc, reaction, velocity_ref, problem->mu, SICONOS_FRICTION_3D_SCALING_CONE_MU,
      SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE, SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX,
      1.0, cblas_dnrm2(n, problem->q, 1));
  double error_ref = err_dyn + err_comp;

  /* Vectorized version (modifies velocity in place). */
  memcpy(velocity, problem->b, m * sizeof(double));
  NM_tgemv(1, problem->H, globalVelocity, 1.0, velocity);
  double error = 0.0;
  int err_info = gfc3d_compute_error(problem, reaction, velocity, globalVelocity, 1e-10,
                                     options, cblas_dnrm2(n, problem->q, 1),
                                     cblas_dnrm2(m, problem->b, 1), &error);

  double rd = rel_diff(error, error_ref);
  printf("gfc3d_compute_error (CONE_MU_PER_CONE_MAX) rel diff = %e (info=%d)\n", rd, err_info);
  if (rd > 1e-12) {
    fprintf(stderr, "gfc3d_compute_error differs from scalar reference!\n");
    info += 1;
  }

  free(reaction);
  free(velocity);
  free(globalVelocity);
  free(velocity_ref);
  free(tmp);
  free(tmp_1);
  free(problem->mu);
  free(problem->q);
  free(problem->b);
  free(problem->M->matrix0);
  problem->M->matrix0 = NULL;
  NM_free(problem->M);
  free(problem->H->matrix0);
  problem->H->matrix0 = NULL;
  NM_free(problem->H);
  free(problem);
  solver_options_delete(options);

  return info;
}

static int test_scaling_invariance(void) {
  int info = 0;
  const int nc = 3;
  const int n = 3 * nc;
  double r[3 * 3];
  double u[3 * 3];
  double mu[3];
  double r_scaled[3 * 3];
  double u_scaled[3 * 3];
  double work[3];
  double tmp[3 * 3];

  srand48(99);
  /* Use moderate magnitudes so that all safeguarded scalings are in the
     homogeneous regime (norms > 1). */
  for (int i = 0; i < n; ++i) {
    r[i] = ((2.0 * drand48()) - 1.0) * 1e2;
    u[i] = ((2.0 * drand48()) - 1.0) * 1e2;
  }
  for (int i = 0; i < nc; ++i) mu[i] = 0.2 + 0.6 * drand48();

  const double rho = 1e4;
  for (int i = 0; i < n; ++i) {
    r_scaled[i] = rho * r[i];
    u_scaled[i] = rho * u[i];
  }

  int norms[] = {SICONOS_FRICTION_3D_SCALING_EUCLIDEAN, SICONOS_FRICTION_3D_SCALING_CONE,
                 SICONOS_FRICTION_3D_SCALING_CONE_MU,
                 SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED};
  const char* norm_names[] = {"EUCLIDEAN", "CONE", "CONE_MU", "WEIGHTED"};
  const int nnorms = sizeof(norms) / sizeof(norms[0]);

  int scopes[] = {SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
                  SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR};
  const char* scope_names[] = {"PER_CONE", "PER_VECTOR"};
  const int nscopes = sizeof(scopes) / sizeof(scopes[0]);

  int reductions[] = {SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX,
                      SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN,
                      SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED};
  const char* reduction_names[] = {"MAX", "MEAN", "MIXED"};
  const int nreductions = sizeof(reductions) / sizeof(reductions[0]);

  for (int in = 0; in < nnorms; ++in) {
    for (int is = 0; is < nscopes; ++is) {
      for (int ir = 0; ir < nreductions; ++ir) {
        double norm_q = cblas_dnrm2(n, r, 1);
        double err_base = fc3d_compute_complementarity_error(
            nc, r, u, mu, norms[in], scopes[is], reductions[ir], 2.0, work, tmp);
        double err_scaled = fc3d_compute_complementarity_error(
            nc, r_scaled, u_scaled, mu, norms[in], scopes[is], reductions[ir], 2.0, work, tmp);
        double rd = rel_diff(err_base, err_scaled);
        char label[64];
        snprintf(label, sizeof(label), "%s_%s_%s", norm_names[in], scope_names[is],
                 reduction_names[ir]);
        printf("scaling invariance %-40s base=%e scaled=%e rel_diff=%e\n", label, err_base,
               err_scaled, rd);
        if (rd > 1e-12) {
          fprintf(stderr, "%s scaling is not invariant under simultaneous scaling!\n", label);
          info += 1;
        }
      }
    }
  }

  return info;
}

int main(void) {
  int info = 0;
  info += test_fc3d_compute_error_scaling();
  info += test_fc3d_compute_error();
  info += test_gfc3d_compute_error();
  info += test_scaling_invariance();
  return info;
}
