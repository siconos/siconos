// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*
 * Drucker-Prager (Mohr-Coulomb) with cohesion: K = { ||sigma_T|| <= eta sigma_N + c }.
 * plasticity_2d_driver translates the problem (sigma~ = sigma + c/eta e_N, q~ = q - M a).
 *
 *  1. M = I, theta = 0: sigma = P_K(-q) = P_{K_eta}(-q + a) - a, cases tension inside K,
 *     apex (pure tension beyond the cohesion) and sliding;
 *  2. footing problem of the Plasticity tests with a cohesion: solution checked with
 *     plasticity_2d_compute_error (cohesion taken into account), compared with the solution
 *     of the translated problem solved by hand, and different from the solution without
 *     cohesion.
 */

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
#include "plasticity_2d_compute_error.h"
#include "plasticity_2d_solvers.h"

static int failures = 0;

static void check(int ok, const char* what, double value) {
  printf("  %-60s %12.4e  %s\n", what, value, ok ? "ok" : "FAILURE");
  if (!ok) failures++;
}

static double dist(int n, const double* x, const double* y) {
  double s = 0.;
  for (int i = 0; i < n; i++) s += (x[i] - y[i]) * (x[i] - y[i]);
  return sqrt(s);
}

static int solve(PlasticityProblem* problem, double* sigma, double* w, double tol,
                 int itermax) {
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = tol;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = itermax;
  int info = plasticity_2d_driver(problem, sigma, w, options);
  printf("  NSGS: info = %d, iterations = %d, error = %.3e\n", info,
         options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  solver_options_delete(options);
  return info;
}

/* ------------------------------------------------------------------------------------------
 */
static void test_identity(void) {
  printf("=== M = I, eta = 0.5, theta = 0, c = 0.5 (a = c/eta = 1) ===\n");
  int nc = 3, n = 9;
  double* M = (double*)calloc(n * n, sizeof(double));
  for (int i = 0; i < n; i++) M[i + i * n] = 1.;
  double* q = (double*)malloc(n * sizeof(double));
  /* -q = (-0.1, 0, 0): tension within the cohesion, sigma = -q
     -q = (-2, 0, 0): beyond the apex, sigma = apex = (-1, 0, 0)
     -q = (1, 2, 0): sliding, sigma~ = P_K((2, 2, 0)) = (2.4, 1.2, 0), sigma = (1.4, 1.2, 0) */
  const double q0[9] = {0.1, 0., 0., 2., 0., 0., -1., -2., 0.};
  const double expected[9] = {-0.1, 0., 0., -1., 0., 0., 1.4, 1.2, 0.};
  memcpy(q, q0, sizeof q0);
  double* eta = (double*)malloc(nc * sizeof(double));
  double* theta = (double*)calloc(nc, sizeof(double));
  double* cohesion = (double*)malloc(nc * sizeof(double));
  for (int c = 0; c < nc; c++) {
    eta[c] = 0.5;
    cohesion[c] = 0.5;
  }
  Plasticity_DruckerPrager_model* model = plasticity_DruckerPrager_model_new(eta, theta);
  model->cohesion = cohesion;
  PlasticityProblem* problem =
      plasticityProblem_new_with_data(3, nc, NM_create_with_data(NM_DENSE, n, n, M), q, model);

  double sigma[9] = {0.}, w[9] = {0.};
  int info = solve(problem, sigma, w, 1e-14, 1000);
  check(info == 0, "converged", 0.);
  check(dist(n, sigma, expected) < 1e-14, "analytic solution", dist(n, sigma, expected));
  check(problem->q == q && memcmp(q, q0, sizeof q0) == 0 && model->cohesion == cohesion,
        "problem restored (q, cohesion)", 0.);
  double error;
  plasticity_2d_compute_error(problem, sigma, w, 1e-14, NULL, 0., &error);
  check(error < 1e-14, "plasticity_2d_compute_error with cohesion", error);
  /* yield condition with cohesion on the sliding cone */
  double f = hypot(sigma[7], sigma[8]) - (0.5 * sigma[6] + 0.5);
  check(fabs(f) < 1e-14, "sliding cone on the yield surface", f);
  plasticityProblem_free(problem);
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 */
static void test_footing(void) {
  printf("=== footing problem (100 cones) with a cohesion ===\n");
  PlasticityProblem* problem =
      plasticity_new_from_filename("./data/plasticity_2d_footing_100.dat");
  if (!problem) {
    check(0, "cannot read the data file", 0.);
    return;
  }
  int nc = problem->numberOfCones, n = 3 * nc;
  double* sigma0 = (double*)calloc(n, sizeof(double));
  double* sigma1 = (double*)calloc(n, sizeof(double));
  double* sigma2 = (double*)calloc(n, sizeof(double));
  double* w = (double*)calloc(n, sizeof(double));

  printf("-- without cohesion\n");
  solve(problem, sigma0, w, 1e-12, 100000);

  /* cohesion of the order of 1% of the mean normal stress */
  double mean_sigmaN = 0.;
  for (int c = 0; c < nc; c++) mean_sigmaN += fabs(sigma0[3 * c]) / nc;
  double* cohesion = (double*)malloc(nc * sizeof(double));
  for (int c = 0; c < nc; c++) cohesion[c] = 0.01 * mean_sigmaN;
  problem->model.drucker_prager->cohesion = cohesion;

  printf("-- with cohesion c = %g\n", cohesion[0]);
  int info = solve(problem, sigma1, w, 1e-12, 100000);
  check(info == 0, "converged", 0.);
  double error;
  plasticity_2d_compute_error(problem, sigma1, w, 1e-12, NULL, 0., &error);
  check(error < 1e-10, "plasticity_2d_compute_error with cohesion", error);
  double d0 = dist(n, sigma1, sigma0) / (1. + cblas_dnrm2(n, sigma0, 1));
  check(d0 > 1e-6, "different from the solution without cohesion", d0);

  /* translated problem solved by hand: q~ = q - M a, then sigma = sigma~ - a */
  printf("-- translated problem without cohesion\n");
  double* a;
  plasticity_DruckerPrager_cohesion_shift(problem, &a);
  PlasticityProblem* translated = plasticity_copy(problem);
  free(translated->model.drucker_prager->cohesion);
  translated->model.drucker_prager->cohesion = NULL;
  NM_gemv(-1.0, translated->M, a, 1.0, translated->q);
  for (int i = 0; i < n; i++) sigma2[i] = a[i];
  solve(translated, sigma2, w, 1e-12, 100000);
  for (int i = 0; i < n; i++) sigma2[i] -= a[i];
  double d = dist(n, sigma1, sigma2) / (1. + cblas_dnrm2(n, sigma2, 1));
  check(d < 1e-10, "same as the translated problem", d);

  free(a);
  plasticityProblem_free(translated);
  plasticityProblem_free(problem);
  free(sigma0);
  free(sigma1);
  free(sigma2);
  free(w);
  printf("\n");
}

/* sigma = 0 is not a solution when q_N >= 0 but eta |q_T| > q_N + theta |q_T| (theta < eta):
 * the trivial-case tests of the driver and of the Newton local solver used q_N >= 0. */
static void test_not_trivial(void) {
  printf("=== q_N >= 0 but sigma = 0 not solution (eta = 0.5, theta = 0) ===\n");
  int ids[] = {0, PLASTICITY_2D_ONECONE_NSN, PLASTICITY_2D_ONECONE_ProjectionOnCone};
  for (int k = 0; k < 3; k++) {
    double* M = (double*)calloc(9, sizeof(double));
    M[0] = M[4] = M[8] = 1.;
    double* q = (double*)malloc(3 * sizeof(double));
    q[0] = 0.1;
    q[1] = 1.;
    q[2] = 0.;
    double* eta = (double*)malloc(sizeof(double));
    double* theta = (double*)calloc(1, sizeof(double));
    eta[0] = 0.5;
    PlasticityProblem* problem =
        plasticityProblem_new_with_data(3, 1, NM_create_with_data(NM_DENSE, 3, 3, M), q,
                                        plasticity_DruckerPrager_model_new(eta, theta));
    double sigma[3] = {0.}, w[3] = {0.};
    SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
    if (ids[k]) solver_options_update_internal(options, 0, ids[k]);
    options->dparam[SICONOS_DPARAM_TOL] = 1e-14;
    int info = plasticity_2d_driver(problem, sigma, w, options);
    /* theta = 0, M = I: sigma = P_K(-q), -q = (-0.1, -1, 0) -> (0.32, -0.16, 0) */
    const double expected[3] = {0.32, -0.16, 0.};
    printf("  internal solver %s\n", ids[k] ? solver_options_id_to_name(ids[k]) : "default");
    check(info == 0 && dist(3, sigma, expected) < 1e-13, "sigma = P_K(-q)",
          dist(3, sigma, expected));
    solver_options_delete(options);
    plasticityProblem_free(problem);
  }
  printf("\n");
}

int main(void) {
  test_not_trivial();
  test_identity();
  test_footing();
  printf("%d check(s) failed\n", failures);
  return failures ? 1 : 0;
}
