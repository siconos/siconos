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
#include "NumericsMatrix.h"
#include "SOC_scaling.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "fc3d_compute_error.h"

static const double MU = 0.5;
static const double DELTA = 1.0e-6;

/* Perturbation directions (non-aligned with the reference states). */
static const double V_R[3] = {1.0, 0.5, -0.3};
static const double V_U[3] = {0.2, -0.7, 0.4};

typedef enum { REGIME_STICKING = 0, REGIME_SEPARATION, REGIME_SLIDING } Regime;

typedef enum { MAGN_BALANCED = 0, MAGN_R_LARGE, MAGN_U_LARGE } Magnitude;

typedef enum { PERT_NONE = 0, PERT_R, PERT_U, PERT_RU } Perturbation;

typedef struct {
  int norm;
  int scope;
  int reduction;
  const char* label;
} ScalingCase;

static const ScalingCase SCALINGS[] = {
    {SICONOS_FRICTION_3D_SCALING_NO, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "no"},
    {SICONOS_FRICTION_3D_SCALING_EUCLIDEAN, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "euc_pc_max"},
    {SICONOS_FRICTION_3D_SCALING_EUCLIDEAN, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "euc_pv_max"},
    {SICONOS_FRICTION_3D_SCALING_CONE, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "cone_pc_max"},
    {SICONOS_FRICTION_3D_SCALING_CONE_MU, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "cmu_pc_max"},
    {SICONOS_FRICTION_3D_SCALING_CONE_MU, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN, "cmu_pv_mean"},
    {SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED, "wgt_pc_mix"},
    {SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR,
     SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "wgt_pv_max"}};
static const int N_SCALINGS = (int)(sizeof(SCALINGS) / sizeof(SCALINGS[0]));

static const char* REGIME_NAMES[] = {"sticking", "separation", "sliding"};
static const char* MAGN_NAMES[] = {"balanced", "r_large", "u_large"};
static const char* PERT_NAMES[] = {"none", "dr", "du", "dr+du"};

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void normalize3(const double v[3], double out[3]) {
  double n = norm3(v);
  if (n < DBL_EPSILON) n = 1.0;
  for (int i = 0; i < 3; ++i) out[i] = v[i] / n;
}

/* Set an exact closed-form solution (r, u) for the chosen regime and magnitude.
 *
 * The problem is constructed with M = I and q = u - r, so that z = r gives
 * w = M z + q = u and the pair (r, u) is an exact solution of the FC3D.
 */
static void set_exact_state(Regime regime, Magnitude magn, double r[3], double u[3]) {
  double R = 1.0, U = 1.0;
  switch (magn) {
    case MAGN_BALANCED:
      R = 1.0;
      U = 1.0;
      break;
    case MAGN_R_LARGE:
      R = 1.0e6;
      U = 1.0e-3;
      break;
    case MAGN_U_LARGE:
      R = 1.0e-3;
      U = 1.0e6;
      break;
  }

  for (int i = 0; i < 3; ++i) {
    r[i] = 0.0;
    u[i] = 0.0;
  }

  switch (regime) {
    case REGIME_STICKING:
      /* r strictly inside the cone, u = 0. */
      r[0] = R;
      break;

    case REGIME_SEPARATION:
      /* r = 0, u_mod in the dual cone (u_N > 0 is enough for μ < 1). */
      u[0] = U;
      break;

    case REGIME_SLIDING:
      /* r on the cone boundary, u_N = 0 and u_T opposite r_T. */
      r[0] = R;
      r[1] = MU * R;
      u[1] = -U;
      break;
  }
}

static void apply_perturbation(const double r[3], const double u[3], Perturbation pert,
                               double rp[3], double up[3]) {
  double vr[3], vu[3];
  normalize3(V_R, vr);
  normalize3(V_U, vu);

  for (int i = 0; i < 3; ++i) {
    rp[i] = r[i];
    up[i] = u[i];
    if (pert == PERT_R || pert == PERT_RU) rp[i] += DELTA * vr[i];
    if (pert == PERT_U || pert == PERT_RU) up[i] += DELTA * vu[i];
  }
}

static FrictionContactProblem* make_identity_problem(const double r[3], const double u[3]) {
  FrictionContactProblem* problem = frictionContactProblem_new();
  problem->numberOfContacts = 1;
  problem->dimension = 3;
  problem->mu = (double*)malloc(sizeof(double));
  problem->mu[0] = MU;
  problem->q = (double*)malloc(3 * sizeof(double));
  for (int i = 0; i < 3; ++i) problem->q[i] = u[i] - r[i];

  NumericsMatrix* M = NM_create_with_data(NM_DENSE, 3, 3, calloc(9, sizeof(double)));
  for (int i = 0; i < 3; ++i) M->matrix0[i * 3 + i] = 1.0;
  problem->M = M;
  return problem;
}

static void free_problem(FrictionContactProblem* problem) {
  if (!problem) return;
  free(problem->mu);
  free(problem->q);
  if (problem->M) {
    free(problem->M->matrix0);
    problem->M->matrix0 = NULL;
    NM_free(problem->M);
  }
  free(problem);
}

/* Compute the per-contact scale for a single cone.
   The scope is ignored because with one contact per-vector and per-cone coincide. */
static double compute_sigma(const double r[3], const double u[3], int norm, int scope,
                            int reduction, double alpha, double norm_q) {
  (void)scope;
  (void)norm_q;
  double u_mod[3];
  double work[1];
  fc3d_compute_modified_local_velocity(1, u, &MU, u_mod, work);

  switch (norm) {
    case SICONOS_FRICTION_3D_SCALING_NO:
      return 1.0;

    case SICONOS_FRICTION_3D_SCALING_EUCLIDEAN:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          return SOC_scaling_euclidean(r, u_mod, 1, 3);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          return SOC_geom_scaling_euclidean(r, u_mod, 1, 3);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          return SOC_mixed_scaling_euclidean(r, u_mod, 1, 3);
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          return SOC_scaling_cone(r, u_mod, 1, 3);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          return SOC_geom_scaling_cone(r, u_mod, 1, 3);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          return SOC_mixed_scaling_cone(r, u_mod, 1, 3);
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE_MU:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          return SOC_scaling_cone_mu(r, u_mod, 1, &MU, 3);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          return SOC_geom_scaling_cone_mu(r, u_mod, 1, &MU, 3);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          return SOC_mixed_scaling_cone_mu(r, u_mod, 1, &MU, 3);
      }
      break;

    case SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED:
      switch (reduction) {
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX:
          return SOC_scaling_cone_weighted(r, u_mod, 1, 3, alpha);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN:
          return SOC_geom_scaling_cone_weighted(r, u_mod, 1, 3, alpha);
        case SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED:
          return SOC_mixed_scaling_cone_weighted(r, u_mod, 1, 3, alpha);
      }
      break;
  }
  return 1.0;
}

static void print_separator(void) {
  printf("+----------------+----------+---------+");
  for (int s = 0; s < N_SCALINGS; ++s) printf("--------------+");
  printf("\n");
}

static void print_header(const char* quantity) {
  print_separator();
  printf("| %-14s | %-8s | %-7s |", "regime", "magn", "pert");
  for (int s = 0; s < N_SCALINGS; ++s) printf(" %12s |", SCALINGS[s].label);
  printf("\n");
  printf("| %-14s | %-8s | %-7s |", "", "", quantity);
  for (int s = 0; s < N_SCALINGS; ++s) printf(" %12s |", "");
  printf("\n");
  print_separator();
}

int main(void) {
  int info = 0;
  SolverOptions options;
  memset(&options, 0, sizeof(options));
  options.iparam = (int*)calloc(OPTIONS_PARAM_SIZE, sizeof(int));
  options.dparam = (double*)calloc(OPTIONS_PARAM_SIZE, sizeof(double));

  /* Table of errors. */
  print_header("error");

  for (int regime_i = 0; regime_i < 3; ++regime_i) {
    for (int magn_i = 0; magn_i < 3; ++magn_i) {
      for (int pert_i = 0; pert_i < 4; ++pert_i) {
        double r[3], u[3], rp[3], up[3];
        set_exact_state((Regime)regime_i, (Magnitude)magn_i, r, u);
        apply_perturbation(r, u, (Perturbation)pert_i, rp, up);

        FrictionContactProblem* problem = make_identity_problem(rp, up);
        double norm_q = cblas_dnrm2(3, problem->q, 1);
        double w[3];

        double errors[N_SCALINGS];
        for (int s = 0; s < N_SCALINGS; ++s) {
          options.iparam[SICONOS_FRICTION_3D_SCALING_NORM] = SCALINGS[s].norm;
          options.iparam[SICONOS_FRICTION_3D_SCALING_SCOPE] = SCALINGS[s].scope;
          options.iparam[SICONOS_FRICTION_3D_SCALING_REDUCTION] = SCALINGS[s].reduction;
          options.dparam[SICONOS_FRICTION_3D_SCALING_ALPHA] = 1.0;
          errors[s] = 0.0;
          fc3d_compute_error(problem, rp, w, /* tolerance */ 1.0, &options, &errors[s]);
        }

        printf("| %-14s | %-8s | %-7s |", REGIME_NAMES[regime_i], MAGN_NAMES[magn_i],
               PERT_NAMES[pert_i]);
        for (int s = 0; s < N_SCALINGS; ++s) printf(" %12.4e |", errors[s]);
        printf("\n");
        fflush(stdout);

        /* Sanity check: unperturbed exact states must give (near) zero error.
           The tolerance is relaxed for large-magnitude states because the
           projection is evaluated in double precision. */
        if (pert_i == PERT_NONE) {
          double scale = fmax(1.0, fmax(cblas_dnrm2(3, rp, 1), cblas_dnrm2(3, up, 1)));
          double tol = 1e-10 * scale;
          for (int s = 0; s < N_SCALINGS; ++s) {
            if (errors[s] > tol) {
              fprintf(stderr, "FAIL: exact state %s/%s/%s has error %e with scaling %s\n",
                      REGIME_NAMES[regime_i], MAGN_NAMES[magn_i], PERT_NAMES[pert_i],
                      errors[s], SCALINGS[s].label);
              info += 1;
            }
          }
        }

        free_problem(problem);
      }
    }
  }
  print_separator();

  /* Table of scales σ. */
  printf("\n");
  print_header("sigma");

  for (int regime_i = 0; regime_i < 3; ++regime_i) {
    for (int magn_i = 0; magn_i < 3; ++magn_i) {
      for (int pert_i = 0; pert_i < 4; ++pert_i) {
        double r[3], u[3], rp[3], up[3];
        set_exact_state((Regime)regime_i, (Magnitude)magn_i, r, u);
        apply_perturbation(r, u, (Perturbation)pert_i, rp, up);

        double q[3];
        for (int i = 0; i < 3; ++i) q[i] = up[i] - rp[i];
        double norm_q = cblas_dnrm2(3, q, 1);

        printf("| %-14s | %-8s | %-7s |", REGIME_NAMES[regime_i], MAGN_NAMES[magn_i],
               PERT_NAMES[pert_i]);
        for (int s = 0; s < N_SCALINGS; ++s) {
          double sigma = compute_sigma(rp, up, SCALINGS[s].norm, SCALINGS[s].scope,
                                       SCALINGS[s].reduction, 1.0, norm_q);
          printf(" %12.4e |", sigma);
        }
        printf("\n");
      }
    }
  }
  print_separator();

  free(options.iparam);
  free(options.dparam);
  return info;
}
