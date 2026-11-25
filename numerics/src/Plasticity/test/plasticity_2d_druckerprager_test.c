/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*
 * Drucker-Prager one-cone solvers (plasticity_2d_druckerprager.h):
 *   sigma in K, -(w + theta ||w_T|| e_N) in N_K(sigma), K = {||sigma_T|| <= eta sigma_N}.
 *
 *  1. return mapping, manufactured solutions: for random non-symmetric M (positive definite
 *     symmetric part), eta, theta in {0, eta/2, eta, 1.5 eta}, a solution of each kind is
 *     built (elastic: sigma inside K, w = 0; apex: sigma = 0, w + theta ||w_T|| e_N in the
 *     dual cone; sliding: sigma = sigma_N (1, eta t), w = lambda (eta - theta, -t)) and
 *     q = w - M sigma. With the manufactured sigma as input (choice among several solutions),
 *     the return mapping must give it back;
 *  2. random local problems: return mapping (residual), semi-smooth Newton without fallback
 *     (convergence, same solution as the return mapping for theta = eta and symmetric M);
 *  3. NSGS on the footing problems with all the Drucker-Prager local solvers: same solution
 *     (tolerance 1e-10: the error is absolute and |sigma| ~ 5e3, 1e-12 is at round-off).
 */

#include "plasticity_2d_druckerprager.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
  printf("  %-66s %12.4e  %s\n", what, value, ok ? "ok" : "FAILURE");
  if (!ok) failures++;
}

static double uniform(void) { return 2. * rand() / RAND_MAX - 1.; }

static double dist3(const double* x, const double* y) {
  return sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
              (x[2] - y[2]) * (x[2] - y[2]));
}

static void random_matrix(double M[9], int symmetric) {
  for (int i = 0; i < 9; i++) M[i] = 0.4 * uniform();
  if (symmetric)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < i; j++) M[i + 3 * j] = M[j + 3 * i];
  M[0] += 2.;
  M[4] += 1. + rand() / (double)RAND_MAX;
  M[8] += 1.5;
}

static double residual(const double M[9], const double q[3], double eta, double theta,
                       const double s[3]) {
  /* rho = 1 / max diag, scale-free enough for these well-conditioned matrices */
  return plasticity_2d_DruckerPrager_natural_map_residual(M, q, eta, theta, 0.5, s);
}

/* ------------------------------------------------------------------------------------------
 */
static void test_manufactured(void) {
  printf("=== return mapping: manufactured elastic / apex / sliding solutions ===\n");
  srand(1);
  const char* kinds[3] = {"elastic", "apex", "sliding"};
  int count[3] = {0, 0, 0}, found[3] = {0, 0, 0};
  double max_dist[3] = {0., 0., 0.};
  for (int k = 0; k < 3000; k++) {
    int kind = k % 3;
    double M[9], q[3], s[3], w[3];
    random_matrix(M, k % 2);
    double eta = 0.2 + 0.8 * rand() / (double)RAND_MAX;
    double theta = eta * 0.5 * ((k / 3) % 4); /* 0, eta/2, eta, 1.5 eta */
    double phi = 2. * M_PI * rand() / (double)RAND_MAX;
    double t[2] = {cos(phi), sin(phi)};
    if (kind == 0) { /* elastic: strictly inside K, w = 0 */
      double sN = 0.5 + rand() / (double)RAND_MAX, r = 0.9 * eta * sN * rand() / RAND_MAX;
      s[0] = sN, s[1] = r * t[0], s[2] = r * t[1];
      w[0] = w[1] = w[2] = 0.;
    } else if (kind == 1) { /* apex: v = w + theta |w_T| e_N strictly in the dual cone */
      s[0] = s[1] = s[2] = 0.;
      double wT = rand() / (double)RAND_MAX;
      /* dual cone: eta |v_T| <= v_N */
      double vN = eta * wT + 0.1 + rand() / (double)RAND_MAX;
      w[0] = vN - theta * wT, w[1] = wT * t[0], w[2] = wT * t[1];
    } else { /* sliding: sigma = sN (1, eta t), w = lambda (eta - theta, -t) */
      double sN = 0.5 + rand() / (double)RAND_MAX, lambda = 0.1 + rand() / (double)RAND_MAX;
      s[0] = sN, s[1] = sN * eta * t[0], s[2] = sN * eta * t[1];
      w[0] = lambda * (eta - theta), w[1] = -lambda * t[0], w[2] = -lambda * t[1];
    }
    for (int i = 0; i < 3; i++)
      q[i] = w[i] - (M[i] * s[0] + M[i + 3] * s[1] + M[i + 6] * s[2]);
    double x[3] = {s[0], s[1], s[2]};
    int status = plasticity_2d_DruckerPrager_return_mapping(M, q, eta, theta, x);
    count[kind]++;
    double d = dist3(x, s) / (1. + sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]));
    if (status == 0 && d < 1e-10) found[kind]++;
    max_dist[kind] = fmax(max_dist[kind], d);
  }
  for (int kind = 0; kind < 3; kind++) {
    char what[128];
    snprintf(what, sizeof what, "%s: %d/%d manufactured solutions recovered", kinds[kind],
             found[kind], count[kind]);
    check(found[kind] == count[kind], what, max_dist[kind]);
  }
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 */
static void test_random_local(void) {
  printf("=== random local problems: return mapping and semi-smooth Newton ===\n");
  srand(2);
  int n = 2000, rm_fail = 0, nsn_fail = 0, compared = 0, plastic = 0;
  double max_rm = 0., max_nsn = 0., max_diff = 0.;
  long it = 0;
  for (int k = 0; k < n; k++) {
    double M[9], q[3];
    int symmetric = k % 2;
    random_matrix(M, symmetric);
    for (int i = 0; i < 3; i++) q[i] = 2. * uniform();
    double eta = 0.2 + 0.8 * rand() / (double)RAND_MAX;
    double theta = eta * 0.5 * ((k / 2) % 3); /* 0, eta/2, eta */
    double s_rm[3] = {0., 0., 0.}, s_nsn[3] = {0., 0., 0.};
    rm_fail += plasticity_2d_DruckerPrager_return_mapping(M, q, eta, theta, s_rm) != 0;
    max_rm = fmax(max_rm, residual(M, q, eta, theta, s_rm));
    if (fabs(hypot(s_rm[1], s_rm[2]) - eta * s_rm[0]) < 1e-12 && s_rm[0] > 1e-12) plastic++;
    Plasticity2DLocalParameters p = {
        .max_iter = 100,
        .tol = 1e-14,
        .rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM,
        .linesearch = PLASTICITY_NSN_LINESEARCH_ARMIJO,
        .linesearch_max_iter = 20,
        .hybrid = 0};
    nsn_fail += plasticity_2d_DruckerPrager_semismooth_newton(M, q, eta, theta, s_nsn, &p);
    it += p.iter_done;
    max_nsn = fmax(max_nsn, residual(M, q, eta, theta, s_nsn));
    if (symmetric && theta == eta) { /* unique solution */
      compared++;
      max_diff = fmax(max_diff, dist3(s_rm, s_nsn));
    }
  }
  printf("  %d problems, %d sliding; Newton: %.2f iterations on average\n", n, plastic,
         (double)it / n);
  check(rm_fail == 0, "return mapping: a solution found", (double)rm_fail);
  check(max_rm < 1e-14, "return mapping: natural-map residual", max_rm);
  check(nsn_fail <= n / 100, "semi-smooth Newton (no fallback): failures <= 1%",
        (double)nsn_fail);
  check(max_diff < 1e-12, "theta = eta, M symmetric: Newton = return mapping", max_diff);
  printf("  (%d failures of the Newton method without fallback, %d unique cases compared)\n",
         nsn_fail, compared);
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 */
static void test_footing(const char* file) {
  printf("=== NSGS on %s ===\n", file);
  PlasticityProblem* problem = plasticity_new_from_filename(file);
  if (!problem) {
    check(0, "cannot read the data file", 0.);
    return;
  }
  int n = 3 * problem->numberOfCones;
  int ids[5] = {0, PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING,
                PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN, PLASTICITY_2D_ONECONE_NSN,
                PLASTICITY_2D_ONECONE_ProjectionOnCone};
  double* ref = (double*)calloc(n, sizeof(double));
  double* s = (double*)calloc(n, sizeof(double));
  double* w = (double*)calloc(n, sizeof(double));
  for (int k = 0; k < 5; k++) {
    SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
    if (ids[k]) solver_options_update_internal(options, 0, ids[k]);
    options->dparam[SICONOS_DPARAM_TOL] = 1e-10; /* absolute error, |sigma| ~ 5e3 */
    options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
    memset(s, 0, n * sizeof(double));
    clock_t t0 = clock();
    int info = plasticity_2d_driver(problem, s, w, options);
    double time = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("  %-52s info %d, %6d iterations, error %.2e, %.3f s\n",
           ids[k] ? solver_options_id_to_name(ids[k]) : "default", info,
           options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU],
           time);
    check(info == 0, "converged", 0.);
    if (k == 0)
      memcpy(ref, s, n * sizeof(double));
    else {
      double d = 0., nr = 0.;
      for (int i = 0; i < n; i++)
        d += (s[i] - ref[i]) * (s[i] - ref[i]), nr += ref[i] * ref[i];
      d = sqrt(d) / (1. + sqrt(nr));
      /* stopped on an absolute error of 1e-10 with M a compliance (~1e-7): the stresses of
       * the different local solvers (iteration paths) differ by ~1e-8 relative */
      check(d < 1e-7, "same solution as the default local solver (relative)", d);
    }
    solver_options_delete(options);
  }
  free(ref);
  free(s);
  free(w);
  plasticityProblem_free(problem);
  printf("\n");
}

int main(void) {
  test_manufactured();
  test_random_local();
  test_footing("./data/plasticity_2d_footing_100.dat");
  test_footing("./data/plasticity_2d_footing_100_theta0.dat");
  test_footing("./data/plasticity_2d_footing_100_theta0.05.dat");
  printf("%d check(s) failed\n", failures);
  return failures ? 1 : 0;
}
