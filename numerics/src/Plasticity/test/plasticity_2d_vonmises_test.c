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
 * Von Mises (J2) plasticity in cone coordinates sigma = (sigma_N, sigma_T1, sigma_T2):
 *   C = { ||sigma_T|| <= R },  R = sqrt(2/3) sigma_y,
 *   w = M sigma + q,  sigma in C,  -w in N_C(sigma)  <=>  sigma = P_C(sigma - rho w).
 *
 * Tests:
 *  1. projection on C;
 *  2. exact one-element solver: analytic cases (M = I), non-symmetric local matrices
 *     (elastic and plastic), 1000 pseudo-random local problems;
 *  3. NSGS solvers (PLASTICITY_2D_NSGS with the default, the exact and the projection local
 *     solvers, PLASTICITY_2D_NSGS_GENERIC) on the Von Mises data files.
 * The solutions are checked with the natural-map residual |sigma - P_C(sigma - rho w)|,
 * rho = 1 / max_i M_ii, w recomputed as M sigma + q.
 */

#include "plasticity_2d_vonmises.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "plasticity_2d_solvers.h"

static int failures = 0;

static void check(int ok, const char* what, double value) {
  printf("  %-60s %12.4e  %s\n", what, value, ok ? "ok" : "FAILURE");
  if (!ok) failures++;
}

/* residual |sigma - P_C(sigma - rho (M sigma + q))| for nc elements, M dense column-major
 * (n x n), rho = 1 / max_i M_ii */
static double vonmises_residual(int nc, const double* M, const double* q,
                                const double* sigma_y, const double* sigma) {
  int n = 3 * nc;
  double dmax = 0.;
  for (int i = 0; i < n; i++) dmax = fmax(dmax, fabs(M[i + i * n]));
  double* w = (double*)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    w[i] = q[i];
    for (int j = 0; j < n; j++) w[i] += M[i + j * n] * sigma[j];
    w[i] /= dmax;
  }
  double err = 0., worktmp[3];
  for (int c = 0; c < nc; c++)
    plasticity_2d_vonmises_unitary_compute_and_add_error(sigma + 3 * c, w + 3 * c, sigma_y[c],
                                                         &err, worktmp);
  free(w);
  return sqrt(err);
}

/* ------------------------------------------------------------------------------------------
 */
/* k-th pseudo-random local problem: M = A A^T + 0.2 I + skew part (positive definite
 * symmetric part), q of amplitude 1e-3 to 1e3 */
static void random_local_problem(int k, double M[9], double q[3], double* sigma_y) {
  double A[9];
  for (int i = 0; i < 9; i++) A[i] = sin(12.9898 * (k + 1) + 78.233 * i);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.;
      for (int l = 0; l < 3; l++) s += A[i + 3 * l] * A[j + 3 * l];
      double skew = 0.3 * sin(3.7 * (k + 1) + 1.3 * i - 2.1 * j) -
                    0.3 * sin(3.7 * (k + 1) + 1.3 * j - 2.1 * i);
      M[i + 3 * j] = s + (i == j ? 0.2 : 0.) + skew;
    }
  double amp = pow(10., (k % 7) - 3);
  for (int i = 0; i < 3; i++) q[i] = amp * sin(4.1 * (k + 1) + 2.7 * i);
  *sigma_y = 0.5 + 0.4 * sin(0.37 * k);
}

static void test_projection(void) {
  printf("=== projection on the Von Mises cylinder ===\n");
  double sigma_y = 2.0, R = plasticity_2d_vonMises_radius(sigma_y);
  check(fabs(R - sqrt(2. / 3.) * 2.0) < 1e-15, "R = sqrt(2/3) sigma_y", R);

  double inside[3] = {-5., 0.3 * R, -0.4 * R};
  double copy[3];
  memcpy(copy, inside, sizeof copy);
  int status = plasticity_2d_projectionOnVonMises(inside, sigma_y);
  check(status == 1 && memcmp(inside, copy, sizeof copy) == 0, "inside: unchanged", 0.);
  check(plasticity_2d_vonMises_check_yield(inside, sigma_y) == 1, "inside: check_yield", 0.);

  double outside[3] = {7., 3., -4.}; /* |sigma_T| = 5 > R */
  status = plasticity_2d_projectionOnVonMises(outside, sigma_y);
  check(status == 0, "outside: projected", 0.);
  check(outside[0] == 7., "outside: sigma_N unchanged", outside[0] - 7.);
  check(fabs(hypot(outside[1], outside[2]) - R) < 1e-15 * R, "outside: |sigma_T| = R",
        hypot(outside[1], outside[2]) - R);
  check(fabs(outside[1] * 4. + outside[2] * 3.) < 1e-15, "outside: radial direction",
        outside[1] * 4. + outside[2] * 3.);
  check(fabs(plasticity_2d_vonMises_equivalent_stress(outside) - sigma_y) < 1e-14,
        "outside: sigma_eq = sigma_y", plasticity_2d_vonMises_equivalent_stress(outside));
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 */
/* Check the optimality conditions of a one-element solution: residual, w_N = 0, and either
 * |sigma_T| < R and w_T = 0, or |sigma_T| = R and w_T = -lambda sigma_T / R, lambda >= 0. */
static int check_local_solution(const double M[9], const double q[3], double sigma_y,
                                const double sigma[3], double tol, int verbose_check) {
  double R = plasticity_2d_vonMises_radius(sigma_y);
  double w[3];
  for (int i = 0; i < 3; i++)
    w[i] = q[i] + M[i] * sigma[0] + M[i + 3] * sigma[1] + M[i + 6] * sigma[2];
  double scale = fmax(1., fmax(fabs(sigma[0]), fmax(R, hypot(q[1], q[2]) + fabs(q[0]))));
  double res = vonmises_residual(1, M, q, &sigma_y, sigma) / scale;
  double nT = hypot(sigma[1], sigma[2]);
  double wN = fabs(w[0]) / scale;
  double compl;
  if (nT < R * (1. - 1e-12)) {
    compl = hypot(w[1], w[2]) / scale; /* elastic: w = 0 */
  } else {
    /* plastic: w_T x sigma_T = 0 and w_T . sigma_T <= 0 */
    double cross = fabs(w[1] * sigma[2] - w[2] * sigma[1]) / (scale * R);
    double dot = (w[1] * sigma[1] + w[2] * sigma[2]) / (scale * R);
    compl = fmax(cross, fmax(dot, 0.)) + fabs(nT - R) / R;
  }
  int ok = (res <= tol) && (wN <= tol) && (compl <= tol);
  if (verbose_check || !ok) {
    printf("    sigma = (% .12e, % .12e, % .12e), |sigma_T| / R = %.15f\n", sigma[0], sigma[1],
           sigma[2], nT / R);
    printf("    residual = %.3e, |w_N| = %.3e, complementarity = %.3e\n", res, wN, compl);
  }
  return ok;
}

static void test_local_solver(void) {
  printf("=== exact one-element solver ===\n");
  double sigma[3];

  /* M = I, q = (-1, 1, 3), sigma_y = 0.5 (data of plasticity_2d_vonmises_example1.dat):
     sigma_N = -q_N = 1, sigma_T = -q_T R / |q_T| */
  const double I3[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  const double q1[3] = {-1., 1., 3.};
  double sigma_y = 0.5, R = plasticity_2d_vonMises_radius(sigma_y);
  int status = plasticity_2d_vonMises_return_mapping(I3, q1, sigma_y, sigma);
  double expected[3] = {1., -R / sqrt(10.), -3. * R / sqrt(10.)};
  double err = fabs(sigma[0] - expected[0]) + fabs(sigma[1] - expected[1]) +
               fabs(sigma[2] - expected[2]);
  check(status == 0 && err < 1e-15, "M = I, plastic: analytic solution", err);
  check(check_local_solution(I3, q1, sigma_y, sigma, 1e-14, 0), "M = I, plastic: conditions",
        0.);

  /* M = I, elastic: sigma = -q */
  const double q2[3] = {2., 0.1, -0.2};
  status = plasticity_2d_vonMises_return_mapping(I3, q2, sigma_y, sigma);
  err = fabs(sigma[0] + q2[0]) + fabs(sigma[1] + q2[1]) + fabs(sigma[2] + q2[2]);
  check(status == 0 && err < 1e-15, "M = I, elastic: sigma = -q", err);

  /* non-symmetric M with coupling between the hydrostatic and deviatoric parts */
  const double M_rows[9] = {2.0, 0.3, -0.2, 0.1, 1.5, 0.4, -0.3, 0.2, 1.2};
  double M[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M[i + 3 * j] = M_rows[3 * i + j];
  const double q3[3] = {-1., 2., -1.5}; /* plastic */
  status = plasticity_2d_vonMises_return_mapping(M, q3, 1.0, sigma);
  check(status == 0 && check_local_solution(M, q3, 1.0, sigma, 1e-13, 1),
        "non-symmetric M, plastic: conditions", 0.);
  const double q4[3] = {-1., 0.05, -0.02}; /* elastic */
  status = plasticity_2d_vonMises_return_mapping(M, q4, 1.0, sigma);
  check(status == 0 && check_local_solution(M, q4, 1.0, sigma, 1e-13, 1),
        "non-symmetric M, elastic: conditions", 0.);

  /* pseudo-random local problems: M = A A^T + 0.2 I + skew part, q of various sizes */
  int nfail = 0, nplastic = 0;
  double worst = 0.;
  for (int k = 0; k < 1000; k++) {
    double Mk[9], qk[3], sy;
    random_local_problem(k, Mk, qk, &sy);
    status = plasticity_2d_vonMises_return_mapping(Mk, qk, sy, sigma);
    if (hypot(sigma[1], sigma[2]) >= plasticity_2d_vonMises_radius(sy) * (1 - 1e-12))
      nplastic++;
    double res = vonmises_residual(1, Mk, qk, &sy, sigma) /
                 fmax(1., fmax(fabs(sigma[0]), hypot(qk[1], qk[2]) + fabs(qk[0])));
    worst = fmax(worst, res);
    if (status == 1 || !check_local_solution(Mk, qk, sy, sigma, 1e-12, 0)) nfail++;
  }
  printf("  (%d plastic cases out of 1000)\n", nplastic);
  check(nfail == 0, "1000 pseudo-random local problems: worst residual", worst);
  printf("\n");
}

/* Iterative local solvers (fixed point, semi-smooth Newton) on the 1000 local problems,
 * compared with the return mapping. */
typedef int (*local_solver_t)(const double*, const double*, double, double*,
                              Plasticity2DVonMisesLocalParameters*);

static void test_iterative_local_solver(const char* name, local_solver_t solver,
                                        Plasticity2DVonMisesLocalParameters p0,
                                        int must_converge) {
  int nconv = 0, iter_sum = 0, nwrong = 0;
  double worst = 0.;
  for (int k = 0; k < 1000; k++) {
    double Mk[9], qk[3], sy, ref[3], sigma[3] = {0., 0., 0.};
    random_local_problem(k, Mk, qk, &sy);
    plasticity_2d_vonMises_return_mapping(Mk, qk, sy, ref);
    Plasticity2DVonMisesLocalParameters p = p0;
    if (solver(Mk, qk, sy, sigma, &p) == 0) {
      nconv++;
      double d = sqrt((sigma[0] - ref[0]) * (sigma[0] - ref[0]) +
                      (sigma[1] - ref[1]) * (sigma[1] - ref[1]) +
                      (sigma[2] - ref[2]) * (sigma[2] - ref[2])) /
                 (1. + sqrt(ref[0] * ref[0] + ref[1] * ref[1] + ref[2] * ref[2]));
      worst = fmax(worst, d);
      if (d > 1e-8) nwrong++;
    }
    iter_sum += p.iter_done;
  }
  printf("  %-44s converged %4d / 1000, mean iterations %7.1f\n", name, nconv,
         iter_sum / 1000.);
  if (must_converge) check(nconv == 1000, "  all converged", 1000 - nconv);
  check(nwrong == 0, "  same solution as the return mapping (worst rel. diff.)", worst);
}

static void test_iterative_local_solvers(void) {
  printf("=== iterative one-element solvers vs return mapping ===\n");
  Plasticity2DVonMisesLocalParameters p = {.max_iter = 1000,
                                           .tol = 1e-13,
                                           .rho_strategy = 0,
                                           .rho = 0.,
                                           .linesearch = 0,
                                           .linesearch_max_iter = 0,
                                           .hybrid = 0};

  p.rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE;
  test_iterative_local_solver("fixed point, adaptive rho", plasticity_2d_vonMises_fixed_point,
                              p, 1);
  p.rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM;
  test_iterative_local_solver("fixed point, rho = 1/|M|", plasticity_2d_vonMises_fixed_point,
                              p, 0);

  p.max_iter = 50;
  p.rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM;
  p.linesearch = PLASTICITY_NSN_LINESEARCH_ARMIJO;
  test_iterative_local_solver("semi-smooth Newton, Armijo",
                              plasticity_2d_vonMises_semismooth_newton, p, 0);
  p.linesearch = PLASTICITY_NSN_LINESEARCH_GOLDSTEINPRICE;
  test_iterative_local_solver("semi-smooth Newton, Goldstein-Price",
                              plasticity_2d_vonMises_semismooth_newton, p, 0);
  p.linesearch = PLASTICITY_NSN_LINESEARCH_NO;
  test_iterative_local_solver("semi-smooth Newton, no line search",
                              plasticity_2d_vonMises_semismooth_newton, p, 0);
  p.linesearch = PLASTICITY_NSN_LINESEARCH_ARMIJO;
  p.rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_CONSTANT;
  p.rho = 1.0;
  test_iterative_local_solver("semi-smooth Newton, Armijo, rho = 1",
                              plasticity_2d_vonMises_semismooth_newton, p, 0);
  p.rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM;
  p.hybrid = 1;
  test_iterative_local_solver("semi-smooth Newton, Armijo, hybrid",
                              plasticity_2d_vonMises_semismooth_newton, p, 1);
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 */
static void solve_file(const char* filename, int solver_id, int internal_solver_id, double tol,
                       int itermax, double* sigma_out, double res_tol) {
  PlasticityProblem* problem = plasticity_new_from_filename(filename);
  if (!problem) {
    check(0, "cannot read the data file", 0.);
    return;
  }
  int nc = problem->numberOfCones, n = 3 * nc;
  double* sigma = (double*)calloc(n, sizeof(double));
  double* w = (double*)calloc(n, sizeof(double));
  SolverOptions* options = solver_options_create(solver_id);
  options->dparam[SICONOS_DPARAM_TOL] = tol;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = itermax;
  if (internal_solver_id) solver_options_update_internal(options, 0, internal_solver_id);
  int info = plasticity_2d_driver(problem, sigma, w, options);

  double* Md = (double*)calloc((size_t)n * n, sizeof(double));
  NumericsMatrix* Mnm = NM_create_with_data(NM_DENSE, n, n, Md);
  NM_to_dense(problem->M, Mnm);
  double res = vonmises_residual(nc, Md, problem->q, problem->model.von_mises->sigma_y, sigma);
  double scale = 1. + cblas_dnrm2(n, sigma, 1);
  printf("  %s, %s / %s: info = %d, iterations = %d, error = %.3e\n", filename,
         solver_options_id_to_name(solver_id),
         internal_solver_id ? solver_options_id_to_name(internal_solver_id) : "default", info,
         options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  check(info == 0, "converged (info = 0)", options->dparam[SICONOS_DPARAM_RESIDU]);
  check(res / scale <= res_tol, "scaled natural-map residual / (1+|sigma|)", res / scale);
  if (sigma_out) memcpy(sigma_out, sigma, n * sizeof(double));

  NM_free(Mnm);
  solver_options_delete(options);
  plasticityProblem_free(problem);
  free(sigma);
  free(w);
}

static void test_nsgs(void) {
  printf("=== NSGS on plasticity_2d_vonmises_example1.dat (M = I, 3 cones) ===\n");
  const char* ex1 = "./data/plasticity_2d_vonmises_example1.dat";
  double sigma[9];
  double R = plasticity_2d_vonMises_radius(0.5);
  double expected[3] = {1., -R / sqrt(10.), -3. * R / sqrt(10.)};

  int ids[][2] = {
      {PLASTICITY_2D_NSGS, 0},
      {PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING},
      {PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_VONMISES_FP},
      {PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_VONMISES_NSN},
      {PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_ProjectionOnCone},
      {PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration},
      {PLASTICITY_2D_NSGS_GENERIC, 0},
      {PLASTICITY_2D_NSGS_GENERIC, PLASTICITY_2D_ONECONE_VONMISES_FP}};
  for (size_t k = 0; k < sizeof ids / sizeof ids[0]; k++) {
    solve_file(ex1, ids[k][0], ids[k][1], 1e-14, 1000, sigma, 1e-14);
    double err = 0.;
    for (int c = 0; c < 3; c++)
      for (int i = 0; i < 3; i++) err = fmax(err, fabs(sigma[3 * c + i] - expected[i]));
    check(err < 1e-13, "analytic solution", err);
  }

  printf("\n=== NSGS on plasticity_2d_vonmises_footing_100.dat (100 cones) ===\n");
  const char* footing = "./data/plasticity_2d_vonmises_footing_100.dat";
  double* s1 = (double*)calloc(300, sizeof(double));
  double* s2 = (double*)calloc(300, sizeof(double));
  solve_file(footing, PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING, 1e-12,
             100000, s1, 1e-10);
  int footing_ids[][2] = {{PLASTICITY_2D_NSGS, 0},
                          {PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_VONMISES_FP},
                          {PLASTICITY_2D_NSGS, PLASTICITY_2D_ONECONE_VONMISES_NSN},
                          {PLASTICITY_2D_NSGS_GENERIC, 0}};
  for (size_t k = 0; k < sizeof footing_ids / sizeof footing_ids[0]; k++) {
    solve_file(footing, footing_ids[k][0], footing_ids[k][1], 1e-12, 100000, s2, 1e-10);
    double d = 0., nrm = 0.;
    for (int i = 0; i < 300; i++) {
      d += (s1[i] - s2[i]) * (s1[i] - s2[i]);
      nrm += s1[i] * s1[i];
    }
    check(sqrt(d) / (1. + sqrt(nrm)) < 1e-8, "same stress as with the return mapping",
          sqrt(d) / (1. + sqrt(nrm)));
  }
  free(s1);
  free(s2);
  printf("\n");
}

int main(void) {
  test_projection();
  test_local_solver();
  test_iterative_local_solvers();
  test_nsgs();
  printf("%d check(s) failed\n", failures);
  return failures ? 1 : 0;
}
