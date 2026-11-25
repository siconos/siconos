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
#include "plasticity_2d_vonmises.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

#include "Plasticity_options.h"
#include "SolverOptions.h"
#include "naming_conventions.h"
#include "numerics_errors.h"
#include "op3x3.h"
#include "solver_registry.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

double plasticity_2d_vonMises_radius(double sigma_y) { return sqrt(2.0 / 3.0) * sigma_y; }

double plasticity_2d_vonMises_equivalent_stress(const double stress[3]) {
  return sqrt(1.5) * hypot(stress[1], stress[2]);
}

int plasticity_2d_vonMises_check_yield(const double stress[3], double sigma_y) {
  return (hypot(stress[1], stress[2]) <= plasticity_2d_vonMises_radius(sigma_y)) ? 1 : 0;
}

int plasticity_2d_projectionOnVonMises(double stress[3], double sigma_y) {
  double R = plasticity_2d_vonMises_radius(sigma_y);
  double nT = hypot(stress[1], stress[2]);
  DEBUG_PRINTF("Von Mises projection: |sigma_T| = %e, R = %e\n", nT, R);
  if (nT <= R) return 1;
  /* radial return of the deviatoric part, the hydrostatic part is unchanged */
  double scale = R / nT;
  stress[1] *= scale;
  stress[2] *= scale;
  return 0;
}

int plasticity_2d_projectionOnVonMisesBall(double stress[3], double sigma_y) {
  double n = sqrt(stress[0] * stress[0] + stress[1] * stress[1] + stress[2] * stress[2]);
  if (n <= sigma_y) return 1;
  double scale = (sigma_y > 0.0) ? sigma_y / n : 0.0;
  for (int i = 0; i < 3; i++) stress[i] *= scale;
  return 0;
}

void plasticity_2d_vonmises_ball_unitary_compute_and_add_error(const double r[3],
                                                               const double u[3],
                                                               double sigma_y, double* error,
                                                               double worktmp[3]) {
  for (int i = 0; i < 3; i++) worktmp[i] = r[i] - u[i];
  plasticity_2d_projectionOnVonMisesBall(worktmp, sigma_y);
  for (int i = 0; i < 3; i++) *error += (r[i] - worktmp[i]) * (r[i] - worktmp[i]);
}

void plasticity_2d_vonmises_unitary_compute_and_add_error(const double r[3], const double u[3],
                                                          double sigma_y, double* error,
                                                          double worktmp[3]) {
  for (int i = 0; i < 3; i++) worktmp[i] = r[i] - u[i];
  plasticity_2d_projectionOnVonMises(worktmp, sigma_y);
  for (int i = 0; i < 3; i++) *error += (r[i] - worktmp[i]) * (r[i] - worktmp[i]);
}

/* x = -(S + mu I)^{-1} b for the 2x2 matrix S. Returns 0 if (S + mu I) is invertible. */
static int vonmises_solve2(const double S[4], double mu, const double b[2], double x[2]) {
  double a = S[0] + mu, d = S[3] + mu; /* S = [S[0] S[2]; S[1] S[3]] (column-major) */
  double det = a * d - S[2] * S[1];
  if (!(fabs(det) > DBL_MIN)) return 1;
  x[0] = -(d * b[0] - S[2] * b[1]) / det;
  x[1] = -(-S[1] * b[0] + a * b[1]) / det;
  return 0;
}

/* x = -(S + mu I)^{-1} b, S (dim x dim, column-major), dim = 2 or 3 */
static int vonmises_shifted_solve(int dim, const double* S, double mu, const double* b,
                                  double* x) {
  if (dim == 2) return vonmises_solve2(S, mu, b, x);
  double A[9];
  for (int i = 0; i < 9; i++) A[i] = S[i];
  A[0] += mu;
  A[4] += mu;
  A[8] += mu;
  for (int i = 0; i < 3; i++) x[i] = -b[i];
  return solve_3x3_gepp(A, x) ? 1 : 0;
}

static double vonmises_normd(int dim, const double* x) {
  double s = 0.0;
  for (int i = 0; i < dim; i++) s += x[i] * x[i];
  return sqrt(s);
}

/* Plastic case of the return mapping on the ball of radius R > 0 of R^dim: x on input is the
 * elastic trial -S^{-1} b with ||x|| > R; on output x(mu) = -(S + mu I)^{-1} b with mu > 0
 * and ||x|| = R. g(mu) = 1/||x(mu)|| - 1/R is increasing (and concave for a symmetric S):
 * Newton from the left of the root, safeguarded by bisection on the bracket [lo, hi].
 * Returns 0, 1 (singular shifted matrix) or 2 (solved to ~1e3 eps only). */
static int vonmises_shifted_norm_equation(int dim, const double* S, const double* b, double R,
                                          double* x) {
  int status = 0;
  double nx = vonmises_normd(dim, x);
  double lo = 0.0, hi = vonmises_normd(dim, b) / R;
  double xt[3];
  int k = 0;
  /* ||x(mu)|| <= ||b|| / mu when the symmetric part of S is positive semi-definite; in
   * general, enlarge hi until ||x(hi)|| <= R */
  while (k++ < 100 &&
         (vonmises_shifted_solve(dim, S, hi, b, xt) || vonmises_normd(dim, xt) > R)) {
    lo = hi;
    hi *= 2.0;
  }
  double mu = lo;
  double f = nx - R;
  for (int it = 0; it < 100; it++) {
    if (vonmises_shifted_solve(dim, S, mu, b, x)) {
      status = 1;
      break;
    }
    nx = vonmises_normd(dim, x);
    f = nx - R;
    DEBUG_PRINTF("  it %d, mu = %.17e, ||x|| - R = %e\n", it, mu, f);
    if (fabs(f) <= 4.0 * DBL_EPSILON * R) break;
    if (f > 0.0)
      lo = mu;
    else
      hi = mu;
    if (hi - lo <= 2.0 * DBL_EPSILON * hi) break;
    /* Newton step on g: dx/dmu = -(S + mu I)^{-1} x */
    double dx[3];
    if (vonmises_shifted_solve(dim, S, mu, x, dx)) {
      mu = 0.5 * (lo + hi);
      continue;
    }
    double dnx = 0.0; /* d||x||/dmu */
    for (int i = 0; i < dim; i++) dnx += x[i] * dx[i];
    dnx /= nx;
    double g = 1.0 / nx - 1.0 / R, dg = -dnx / (nx * nx);
    double mu_new = (dg != 0.0) ? mu - g / dg : 0.5 * (lo + hi);
    if (!(mu_new > lo && mu_new < hi)) mu_new = 0.5 * (lo + hi);
    mu = mu_new;
  }
  if (status == 0 && fabs(f) > 1e3 * DBL_EPSILON * R) status = 2;
  /* on the yield surface */
  nx = vonmises_normd(dim, x);
  if (nx > 0.0)
    for (int i = 0; i < dim; i++) x[i] *= R / nx;
  return status;
}

int plasticity_2d_vonMises_return_mapping(const double M[9], const double q[3], double sigma_y,
                                          double stress[3]) {
  double m = M[0];
  if (!(m > 0.0)) return 1;
  /* M = [m, M_NT; M_TN, M_TT] (column-major) */
  const double mNT[2] = {M[3], M[6]};
  const double mTN[2] = {M[1], M[2]};
  /* Schur complement S = M_TT - M_TN M_NT / m and q~ = q_T - M_TN q_N / m */
  const double S[4] = {M[4] - mTN[0] * mNT[0] / m, M[5] - mTN[1] * mNT[0] / m,
                       M[7] - mTN[0] * mNT[1] / m, M[8] - mTN[1] * mNT[1] / m};
  const double qt[2] = {q[1] - mTN[0] * q[0] / m, q[2] - mTN[1] * q[0] / m};
  double R = plasticity_2d_vonMises_radius(sigma_y);
  double x[2];
  int status = 0;

  /* elastic trial: S x = -q~ */
  if (vonmises_solve2(S, 0.0, qt, x)) return 1;
  double nx = hypot(x[0], x[1]);
  DEBUG_PRINTF("Von Mises local solve: elastic trial |sigma_T| = %e, R = %e\n", nx, R);

  if (nx > R) {
    if (R <= 0.0) {
      x[0] = x[1] = 0.0;
    } else {
      status = vonmises_shifted_norm_equation(2, S, qt, R, x);
    }
  }
  stress[1] = x[0];
  stress[2] = x[1];
  /* w_N = 0 */
  stress[0] = -(q[0] + mNT[0] * x[0] + mNT[1] * x[1]) / m;
  return status;
}

int plasticity_2d_vonMises_ball_return_mapping(const double M[9], const double q[3],
                                               double sigma_y, double stress[3]) {
  double x[3];
  /* elastic trial: M x = -q */
  if (vonmises_shifted_solve(3, M, 0.0, q, x)) return 1;
  int status = 0;
  if (vonmises_normd(3, x) > sigma_y) {
    if (sigma_y <= 0.0)
      x[0] = x[1] = x[2] = 0.0;
    else
      status = vonmises_shifted_norm_equation(3, M, q, sigma_y, x);
  }
  for (int i = 0; i < 3; i++) stress[i] = x[i];
  return status;
}

/* ---------------------------------------------------------------------------------------- */
/* Iterative solvers on the natural map                                                      */
/* ---------------------------------------------------------------------------------------- */

static void vonmises_w(const double M[9], const double q[3], const double s[3], double w[3]) {
  for (int i = 0; i < 3; i++) w[i] = q[i] + M[i] * s[0] + M[i + 3] * s[1] + M[i + 6] * s[2];
}

static double vonmises_norm3(const double a[3]) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/* default rho = 1 / ||M||_2 */
static double vonmises_default_rho(const double M[9]) {
  double n2 = plasticity_2d_local_spectral_norm(M);
  return (n2 > 0.0) ? 1.0 / n2 : 1.0;
}

/* projection on the cylinder (plane strain) or on the ball */
static void vonmises_project(int ball, double s[3], double sigma_y) {
  if (ball)
    plasticity_2d_projectionOnVonMisesBall(s, sigma_y);
  else
    plasticity_2d_projectionOnVonMises(s, sigma_y);
}

/* Phi = s - P_C(s - rho w); returns |Phi|, y = s - rho w (before projection) */
static double vonmises_natural_map(int ball, const double M[9], const double q[3],
                                   double sigma_y, double rho, const double s[3],
                                   double Phi[3], double y[3]) {
  double w[3];
  vonmises_w(M, q, s, w);
  for (int i = 0; i < 3; i++) Phi[i] = y[i] = s[i] - rho * w[i];
  vonmises_project(ball, Phi, sigma_y);
  for (int i = 0; i < 3; i++) Phi[i] = s[i] - Phi[i];
  return vonmises_norm3(Phi);
}

double plasticity_2d_vonMises_natural_map_residual(const double M[9], const double q[3],
                                                   double sigma_y, double rho,
                                                   const double stress[3]) {
  double Phi[3], y[3];
  return vonmises_natural_map(0, M, q, sigma_y, rho, stress, Phi, y);
}

static double vonmises_residual(int ball, const double M[9], const double q[3], double sigma_y,
                                double rho, const double stress[3]) {
  double Phi[3], y[3];
  return vonmises_natural_map(ball, M, q, sigma_y, rho, stress, Phi, y);
}

int plasticity_2d_vonMises_fixed_point(const double M[9], const double q[3], double sigma_y,
                                       double stress[3],
                                       Plasticity2DVonMisesLocalParameters* p) {
  const double tau = 2.0 / 3.0, tauinv = 3.0 / 2.0, L = 0.9, Lmin = 0.3;
  double rho0 = vonmises_default_rho(M);
  double rho = rho0;
  if (p->rho_strategy == PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_CONSTANT && p->rho > 0.0)
    rho = p->rho;
  int adaptive = (p->rho_strategy == PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE);
  int ls_max = (p->linesearch_max_iter > 0) ? p->linesearch_max_iter : 10;

  double res = vonmises_residual(p->ball, M, q, sigma_y, rho0, stress);
  int it = 0;
  while (res > p->tol * (1.0 + vonmises_norm3(stress)) && it < p->max_iter) {
    it++;
    double w[3], s_new[3];
    vonmises_w(M, q, stress, w);
    if (adaptive) {
      /* decrease rho until rho ||w(s+) - w(s)|| <= L ||s+ - s|| */
      double rho_k = rho / tau, a1 = 0.0, a2 = 0.0;
      int success = 0;
      for (int ls = 0; ls < ls_max && !success; ls++) {
        rho_k *= tau;
        double w_new[3], dw[3], ds[3];
        for (int i = 0; i < 3; i++) s_new[i] = stress[i] - rho_k * w[i];
        vonmises_project(p->ball, s_new, sigma_y);
        vonmises_w(M, q, s_new, w_new);
        for (int i = 0; i < 3; i++) {
          dw[i] = w_new[i] - w[i];
          ds[i] = s_new[i] - stress[i];
        }
        a1 = vonmises_norm3(dw);
        a2 = vonmises_norm3(ds);
        success = (rho_k * a1 <= L * a2);
      }
      rho = (rho_k * a1 < Lmin * a2) ? rho_k * tauinv : rho_k;
    } else {
      for (int i = 0; i < 3; i++) s_new[i] = stress[i] - rho * w[i];
      vonmises_project(p->ball, s_new, sigma_y);
    }
    for (int i = 0; i < 3; i++) stress[i] = s_new[i];
    res = vonmises_residual(p->ball, M, q, sigma_y, rho0, stress);
    DEBUG_PRINTF("Von Mises FP: it %d, rho = %e, residual = %e\n", it, rho, res);
  }
  p->iter_done = it;
  p->residual = res;
  return (res > p->tol * (1.0 + vonmises_norm3(stress))) ? 1 : 0;
}

/* generalized Jacobian of the projection on C at y (column-major) */
static void vonmises_projection_jacobian(int ball, const double y[3], double sigma_y,
                                         double J[9]) {
  if (ball) {
    double n = vonmises_norm3(y);
    for (int i = 0; i < 9; i++) J[i] = 0.0;
    if (n <= sigma_y) {
      J[0] = J[4] = J[8] = 1.0;
    } else {
      /* (r/||y||) (I - u u^T), u = y/||y|| */
      double f = sigma_y / n;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          J[i + 3 * j] = f * ((i == j ? 1.0 : 0.0) - y[i] * y[j] / (n * n));
    }
    return;
  }
  double R = plasticity_2d_vonMises_radius(sigma_y);
  double nT = hypot(y[1], y[2]);
  for (int i = 0; i < 9; i++) J[i] = 0.0;
  J[0] = 1.0;
  if (nT <= R) {
    J[4] = J[8] = 1.0;
  } else {
    double f = R / nT, n1 = y[1] / nT, n2 = y[2] / nT;
    J[4] = f * (1.0 - n1 * n1);
    J[8] = f * (1.0 - n2 * n2);
    J[5] = J[7] = -f * n1 * n2;
  }
}

/* natural map and generalized Jacobian of the Von Mises problem (shared Newton) */
typedef struct {
  const double* M;
  const double* q;
  double sigma_y;
  int ball;
} VonMisesData;

static double vonmises_newton_map(void* data, double rho, const double s[3], double Phi[3],
                                  double H[9]) {
  VonMisesData* d = (VonMisesData*)data;
  double y[3];
  double nPhi = vonmises_natural_map(d->ball, d->M, d->q, d->sigma_y, rho, s, Phi, y);
  if (H) {
    /* H = I - J_P(y) (I - rho M) */
    double J[9], A[9];
    vonmises_projection_jacobian(d->ball, y, d->sigma_y, J);
    for (int i = 0; i < 9; i++) A[i] = -rho * d->M[i];
    A[0] += 1.0;
    A[4] += 1.0;
    A[8] += 1.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double sum = 0.0;
        for (int k = 0; k < 3; k++) sum += J[i + 3 * k] * A[k + 3 * j];
        H[i + 3 * j] = (i == j ? 1.0 : 0.0) - sum;
      }
  }
  return nPhi;
}

static int vonmises_newton_exact(void* data, double s[3]) {
  VonMisesData* d = (VonMisesData*)data;
  int status = d->ball ? plasticity_2d_vonMises_ball_return_mapping(d->M, d->q, d->sigma_y, s)
                       : plasticity_2d_vonMises_return_mapping(d->M, d->q, d->sigma_y, s);
  return (status == 1) ? 1 : 0;
}

int plasticity_2d_vonMises_semismooth_newton(const double M[9], const double q[3],
                                             double sigma_y, double stress[3],
                                             Plasticity2DVonMisesLocalParameters* p) {
  VonMisesData data = {M, q, sigma_y, p->ball};
  return plasticity_2d_local_semismooth_newton(vonmises_newton_map, vonmises_newton_exact,
                                               &data, M, stress, p);
}

/* ===========================================================================
 * Von Mises (J2) local solvers, cone coordinates (SolverOptions interface)
 * =========================================================================== */

int plasticity_2d_projectionOnVonMises_solve(PlasticityProblem* localproblem,
                                             double* stress_local, SolverOptions* options) {
  (void)options;
  /* one projection step: sigma <- P_C(sigma - rho (M sigma + q)), rho = 1 / M_NN */
  double* MLocal = localproblem->M->dense;
  double* qLocal = localproblem->q;
  double sigma_y = localproblem->model.von_mises->sigma_y[0];
  double rho = 1.0 / MLocal[0];
  int ball = (localproblem->model.von_mises->hypothesis == PLASTICITY_VONMISES_BALL);
  double w[3];
  for (int i = 0; i < 3; i++)
    w[i] = qLocal[i] + MLocal[i] * stress_local[0] + MLocal[i + 3] * stress_local[1] +
           MLocal[i + 6] * stress_local[2];
  for (int i = 0; i < 3; i++) stress_local[i] -= rho * w[i];
  vonmises_project(ball, stress_local, sigma_y);
  return 0;
}

int plasticity_2d_vonMises_return_mapping_solve(PlasticityProblem* localproblem,
                                                double* stress_local, SolverOptions* options) {
  double sigma_y = localproblem->model.von_mises->sigma_y[0];
  int status = (localproblem->model.von_mises->hypothesis == PLASTICITY_VONMISES_BALL)
                   ? plasticity_2d_vonMises_ball_return_mapping(
                         localproblem->M->dense, localproblem->q, sigma_y, stress_local)
                   : plasticity_2d_vonMises_return_mapping(
                         localproblem->M->dense, localproblem->q, sigma_y, stress_local);
  if (options) SET_SOLVER_ITER_DONE(options, 1);
  /* status 2: scalar equation solved to ~1e3 eps only, the solution is still usable */
  return (status == 1) ? 1 : 0;
}

/* parameters of the iterative Von Mises local solvers from the local solver options */
static void vonmises_parameters(PlasticityProblem* localproblem, SolverOptions* options,
                                Plasticity2DVonMisesLocalParameters* p) {
  p->ball = (localproblem->model.von_mises->hypothesis == PLASTICITY_VONMISES_BALL);
  p->max_iter = SOLVER_MAX_ITER(options);
  p->tol = SOLVER_TOL(options);
  p->rho_strategy = options->iparam[PLASTICITY_NSN_RHO_STRATEGY];
  p->rho = options->dparam[PLASTICITY_NSN_RHO];
  p->linesearch = options->iparam[PLASTICITY_NSN_LINESEARCH];
  p->linesearch_max_iter = options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER];
  p->hybrid =
      (options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] != PLASTICITY_NSN_HYBRID_STRATEGY_NO);
  p->iter_done = 0;
  p->residual = 0.0;
}

int plasticity_2d_vonMises_fixed_point_solve(PlasticityProblem* localproblem,
                                             double* stress_local, SolverOptions* options) {
  Plasticity2DVonMisesLocalParameters p;
  vonmises_parameters(localproblem, options, &p);
  int info = plasticity_2d_vonMises_fixed_point(localproblem->M->dense, localproblem->q,
                                                localproblem->model.von_mises->sigma_y[0],
                                                stress_local, &p);
  SET_SOLVER_ITER_DONE(options, p.iter_done);
  SET_SOLVER_RESIDUAL(options, p.residual);
  return info;
}

int plasticity_2d_vonMises_nsn_solve(PlasticityProblem* localproblem, double* stress_local,
                                     SolverOptions* options) {
  Plasticity2DVonMisesLocalParameters p;
  vonmises_parameters(localproblem, options, &p);
  int info = plasticity_2d_vonMises_semismooth_newton(
      localproblem->M->dense, localproblem->q, localproblem->model.von_mises->sigma_y[0],
      stress_local, &p);
  SET_SOLVER_ITER_DONE(options, p.iter_done);
  SET_SOLVER_RESIDUAL(options, p.residual);
  return info;
}

void plasticity_2d_vonmises_set_default(SolverOptions* options) { (void)options; }

void plasticity_2d_vonmises_fp_set_default(SolverOptions* options) {
  options->iparam[PLASTICITY_NSN_RHO_STRATEGY] =
      PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE;
  options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER] = 10; /* rho adaptation steps */
  options->dparam[PLASTICITY_NSN_RHO] = 0.0;                /* 0: 1 / ||M||_2 */
}

void plasticity_2d_vonmises_nsn_set_default(SolverOptions* options) {
  options->iparam[PLASTICITY_NSN_RHO_STRATEGY] =
      PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM;
  options->iparam[PLASTICITY_NSN_LINESEARCH] = PLASTICITY_NSN_LINESEARCH_ARMIJO;
  options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER] = 20;
  /* fall back to the return mapping if the Newton method fails */
  options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] =
      PLASTICITY_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP;
  options->dparam[PLASTICITY_NSN_RHO] = 0.0;
}

/* Registry wrappers of the Von Mises local solvers */
static int plasticity_2d_vonmises_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_vonmises_return_mapping_wrap(void* localproblem, double* reaction,
                                                      double* velocity,
                                                      SolverOptions* options) {
  (void)velocity;
  return plasticity_2d_vonMises_return_mapping_solve((PlasticityProblem*)localproblem,
                                                     reaction, options);
}

static int plasticity_2d_vonmises_fp_wrap(void* localproblem, double* reaction,
                                          double* velocity, SolverOptions* options) {
  (void)velocity;
  return plasticity_2d_vonMises_fixed_point_solve((PlasticityProblem*)localproblem, reaction,
                                                  options);
}

static int plasticity_2d_vonmises_nsn_wrap(void* localproblem, double* reaction,
                                           double* velocity, SolverOptions* options) {
  (void)velocity;
  return plasticity_2d_vonMises_nsn_solve((PlasticityProblem*)localproblem, reaction, options);
}

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING,
                "PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING",
                "Von Mises (J2), one element: return mapping (exact solution)",
                plasticity_2d_vonmises_init_wrap, plasticity_2d_vonmises_return_mapping_wrap,
                NULL, NULL, plasticity_2d_vonmises_set_default, 1, /* default_max_iter */
                1e-14,                                             /* default_tol */
                1);                                                /* is_local */

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_VONMISES_FP, "PLASTICITY_2D_ONECONE_VONMISES_FP",
                "Von Mises (J2), one element: fixed point on the natural map",
                plasticity_2d_vonmises_init_wrap, plasticity_2d_vonmises_fp_wrap, NULL, NULL,
                plasticity_2d_vonmises_fp_set_default, 1000, /* default_max_iter */
                1e-14,                                       /* default_tol */
                1);                                          /* is_local */

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_VONMISES_NSN, "PLASTICITY_2D_ONECONE_VONMISES_NSN",
                "Von Mises (J2), one element: semi-smooth Newton on the natural map",
                plasticity_2d_vonmises_init_wrap, plasticity_2d_vonmises_nsn_wrap, NULL, NULL,
                plasticity_2d_vonmises_nsn_set_default, 50, /* default_max_iter */
                1e-14,                                      /* default_tol */
                1);                                         /* is_local */
