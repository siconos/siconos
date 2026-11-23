/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

/*
 * Helpers shared by the MixedProblem tests: build and solve a mixed problem from a dense
 * column-major matrix, and check the solution independently of the MixedProblems code
 * (natural-map residual of each element computed with the unitary error functions of the
 * single-model solvers, w recomputed as M z + q).
 */

#ifndef MIXED_TEST_TOOLS_H
#define MIXED_TEST_TOOLS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MixedProblem.h"
#include "MixedProblem_NSGS.h"
#include "MixedProblems_options.h"
#include "NumericsMatrix.h"
#include "SolverOptions.h"
#include "fc2d_compute_error.h"
#include "fc3d_compute_error.h"
#include "plasticity_2d_compute_error.h"
#include "plasticity_2d_vonmises_projection.h"
#include "rolling_friction_2d_compute_error.h"
#include "rolling_friction_3d_compute_error.h"

/* tolerance of the mixed NSGS (full error, relative to |q|) */
#define MIXED_TOL 1e-10
#define MIXED_MAX_ITER 20000
/* tolerance of the independent checks */
#define CHECK_TOL 1e-8

/* Column-major n x n copy of a row-major array. */
static inline double* colmajor_from_rowmajor(int n, const double* A_rows) {
  double* A = (double*)malloc((size_t)n * n * sizeof(double));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) A[i + j * n] = A_rows[i * n + j];
  return A;
}

/* w = M z + q, M dense column-major. */
static inline void compute_w(int n, const double* M, const double* q, const double* z,
                             double* w) {
  for (int i = 0; i < n; i++) {
    w[i] = q[i];
    for (int j = 0; j < n; j++) w[i] += M[i + j * n] * z[j];
  }
}

static inline double norm2(int n, const double* x) {
  double s = 0.;
  for (int i = 0; i < n; i++) s += x[i] * x[i];
  return sqrt(s);
}

static inline double dist2(int n, const double* x, const double* y) {
  double s = 0.;
  for (int i = 0; i < n; i++) s += (x[i] - y[i]) * (x[i] - y[i]);
  return sqrt(s);
}

/* Description of a block: model, first row, dimension, number of elements and the
   per-element parameters:
     friction 2D/3D: p1 = mu;  rolling friction: p1 = mu, p2 = mu_r;
     Mohr-Coulomb (Drucker-Prager): p1 = eta, p2 = theta;  Von Mises: p1 = sigma_y. */
typedef struct {
  MixedModelType model;
  int start;
  int dim;
  int ne;
  double* p1;
  double* p2;
} BlockSpec;

/* Natural-map residual |z - P_K(z - rho (w + De Saxce/dilatancy term))| of the whole
   problem, with w = M z + q recomputed here and rho = 1 / max_i M_ii, so that z and rho w
   have the same scale (cones and De Saxce terms are positively homogeneous; the Von Mises
   term uses the same rho). */
static inline double mixed_residual(int n, const double* M, const double* q, int nb,
                                    const BlockSpec* blocks, const double* z) {
  double* w = (double*)malloc(n * sizeof(double));
  compute_w(n, M, q, z, w);
  double dmax = 0.;
  for (int i = 0; i < n; i++) dmax = fmax(dmax, fabs(M[i + i * n]));
  for (int i = 0; i < n; i++) w[i] /= dmax;
  double err = 0., worktmp[5];
  for (int b = 0; b < nb; b++) {
    const BlockSpec* B = &blocks[b];
    for (int e = 0; e < B->ne; e++) {
      int idx = B->start + e * B->dim;
      double* ze = (double*)z + idx;
      double* we = w + idx;
      switch (B->model) {
        case MIXED_MODEL_FRICTION_2D:
          fc2d_unitary_compute_and_add_error(ze, we, B->p1[e], &err, worktmp);
          break;
        case MIXED_MODEL_FRICTION_3D:
          fc3d_unitary_compute_and_add_error(ze, we, B->p1[e], &err, worktmp);
          break;
        case MIXED_MODEL_ROLLING_FRICTION_2D:
          rolling_friction_2d_unitary_compute_and_add_error(ze, we, B->p1[e], B->p2[e], &err,
                                                            worktmp);
          break;
        case MIXED_MODEL_ROLLING_FRICTION_3D:
          rolling_friction_3d_unitary_compute_and_add_error(ze, we, B->p1[e], B->p2[e], &err,
                                                            worktmp);
          break;
        case MIXED_MODEL_DRUCKER_PRAGER:
          plasticity_2d_unitary_compute_and_add_error(ze, we, B->p1[e], B->p2[e], &err,
                                                      worktmp);
          break;
        case MIXED_MODEL_VON_MISES:
          for (int i = 0; i < 3; i++) worktmp[i] = ze[i] - we[i];
          plasticity_2d_projectionOnVonMises(worktmp, B->p1[e]);
          for (int i = 0; i < 3; i++) err += (ze[i] - worktmp[i]) * (ze[i] - worktmp[i]);
          break;
        default:
          err = NAN;
          break;
      }
    }
  }
  free(w);
  return sqrt(err);
}

/* Solve with the mixed NSGS; M dense column-major (not modified), given to the solver with
   the storage `storage` (NM_DENSE or NM_SPARSE). Returns the info of the driver. */
static inline int solve_mixed(int n, double* M, double* q, int nb, const BlockSpec* blocks,
                              double* z, double* w, NM_types storage) {
  NumericsMatrix* Md = NM_create_with_data(NM_DENSE, n, n, M); /* M is not owned */
  NumericsMatrix* NM = (storage == NM_SPARSE) ? NM_copy_to_sparse(Md, 0.) : Md;
  MixedProblem* problem = mixedProblem_new_with_data(nb, NM, q);
  for (int b = 0; b < nb; b++) {
    const BlockSpec* B = &blocks[b];
    mixedProblem_set_block(problem, b, B->start, B->dim, B->ne, B->model);
    switch (B->model) {
      case MIXED_MODEL_FRICTION_2D:
        mixedProblem_set_friction2d_params(problem, b, B->p1);
        break;
      case MIXED_MODEL_FRICTION_3D:
        mixedProblem_set_friction3d_params(problem, b, B->p1);
        break;
      case MIXED_MODEL_ROLLING_FRICTION_2D:
        mixedProblem_set_rolling_friction2d_params(problem, b, B->p1, B->p2);
        break;
      case MIXED_MODEL_ROLLING_FRICTION_3D:
        mixedProblem_set_rolling_friction3d_params(problem, b, B->p1, B->p2);
        break;
      case MIXED_MODEL_DRUCKER_PRAGER:
        mixedProblem_set_druckerprager_params(problem, b, B->p1, B->p2);
        break;
      case MIXED_MODEL_VON_MISES:
        mixedProblem_set_vonmises_params(problem, b, B->p1);
        break;
      default:
        break;
    }
  }

  SolverOptions* options = (SolverOptions*)malloc(sizeof(SolverOptions));
  mixedProblem_NSGS_set_default_solver_options(problem, options);
  options->iparam[SICONOS_MIXED_IPARAM_MAX_ITER] = MIXED_MAX_ITER;
  options->dparam[SICONOS_MIXED_DPARAM_TOL] = MIXED_TOL;

  int info = mixedProblem_driver(problem, z, w, options);
  printf("  mixed NSGS (%s M): info = %d, iterations = %d, error = %.3e\n",
         storage == NM_SPARSE ? "sparse" : "dense", info,
         options->iparam[SICONOS_MIXED_IPARAM_ITER_DONE],
         options->dparam[SICONOS_MIXED_DPARAM_RESIDU]);

  solver_options_delete(options);
  mixedProblem_free(problem);
  if (NM != Md) NM_free(NM);
  Md->dense = NULL;
  NM_free(Md);
  return info;
}

/* Check the solver info, the residual, w and (optionally) the expected solution.
   Returns 0 on success. */
static inline int check_solution(const char* name, int solver_info, int n, const double* M,
                                 const double* q, int nb, const BlockSpec* blocks,
                                 const double* z, const double* w, const double* z_expected,
                                 double z_tol) {
  double* w_true = (double*)malloc(n * sizeof(double));
  compute_w(n, M, q, z, w_true);
  double res = mixed_residual(n, M, q, nb, blocks, z);
  double scale = 1. + norm2(n, z);
  double w_err = dist2(n, w, w_true) / (1. + norm2(n, w_true));
  int ok = 1;

  printf("  residual = %.3e (/(1+|z|): %.3e), |w - (Mz+q)| / (1+|Mz+q|) = %.3e\n", res,
         res / scale, w_err);
  if (solver_info != 0) {
    printf("  FAILURE: the solver did not converge (info = %d)\n", solver_info);
    ok = 0;
  }
  /* written as !(x <= tol) so that NaN fails */
  if (!(res / scale <= CHECK_TOL)) {
    printf("  FAILURE: natural-map residual too large\n");
    ok = 0;
  }
  if (!(w_err <= CHECK_TOL)) {
    printf("  FAILURE: returned w differs from M z + q\n");
    ok = 0;
  }
  if (z_expected) {
    double z_err = dist2(n, z, z_expected) / (1. + norm2(n, z_expected));
    printf("  |z - z_expected| / (1+|z_expected|) = %.3e\n", z_err);
    if (!(z_err <= z_tol)) {
      printf("  FAILURE: z differs from the expected solution\n");
      ok = 0;
    }
  }
  if (!ok || n <= 9) {
    int nprint = n < 12 ? n : 12;
    for (int i = 0; i < nprint; i++)
      printf("    z[%2d] = % .10e   w[%2d] = % .10e   Mz+q = % .10e\n", i, z[i], i, w[i],
             w_true[i]);
    if (nprint < n) printf("    ...\n");
  }
  printf("Test %s: %s\n\n", name, ok ? "PASSED" : "FAILED");
  free(w_true);
  return ok ? 0 : 1;
}

#endif /* MIXED_TEST_TOOLS_H */
