// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*
 * Scaling of the error of plasticity_2d_compute_error (as in fc3d_compute_error):
 *  1. no scaling = sqrt(sum of the unitary natural-map errors) (Drucker-Prager with
 *     dilatancy and cohesion, Von Mises);
 *  2. default (Euclidean, per cone, MAX: s_e = min(max(1, |r_e|), max(1, |v_e|))) and MIXED
 *     (s_e = max(1, |r_e| + |v_e|)), PER_VECTOR, computed again here;
 *  3. homogeneity: with (M, lambda q) and lambda sigma, the absolute error is multiplied by
 *     lambda, the MIXED error is unchanged;
 *  4. rho_e = 1/||M_ee||: the error does not depend on the units of the strain rate
 *     ((kappa M, kappa q)); with rho = 1 it does;
 *  5. footing problem (|sigma| ~ 5e3, |w| < 4e-4): with rho_e and the MIXED scaling the NSGS
 *     reaches a tolerance of 1e-13 (the absolute error stalls at ~2e-12, round-off).
 */

#include "plasticity_2d_compute_error.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "SolverOptions.h"
#include "plasticity_2d_solvers.h"
#include "plasticity_2d_vonmises.h"

#define NC 5
#define N (3 * NC)

static int failures = 0;

static void check(int ok, const char* what, double value) {
  printf("  %-66s %12.4e  %s\n", what, value, ok ? "ok" : "FAILURE");
  if (!ok) failures++;
}

static double uniform(void) { return 2. * rand() / RAND_MAX - 1.; }
static double norm3(const double* x) { return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]); }

static PlasticityProblem* random_problem(int drucker_prager, double lambda) {
  srand(7);
  double* M = (double*)calloc(N * N, sizeof(double));
  double* q = (double*)malloc(N * sizeof(double));
  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      if (abs(i / 3 - j / 3) <= 1) M[i + j * N] = (i / 3 == j / 3 ? 0.3 : 0.1) * uniform();
  for (int i = 0; i < N; i++) M[i + i * N] += 2.;
  for (int i = 0; i < N; i++) q[i] = lambda * 3. * uniform();
  double* p1 = (double*)malloc(NC * sizeof(double));
  double* p2 = (double*)malloc(NC * sizeof(double));
  double* p3 = (double*)malloc(NC * sizeof(double));
  for (int e = 0; e < NC; e++) {
    p1[e] = drucker_prager ? 0.3 + 0.1 * e : lambda * (0.5 + 0.2 * e);
    p2[e] = 0.1 * (e % 2);
    p3[e] = lambda * 0.05 * e;
  }
  PlasticityProblem* problem;
  if (drucker_prager) {
    Plasticity_DruckerPrager_model* model = plasticity_DruckerPrager_model_new(p1, p2);
    model->cohesion = p3;
    problem = plasticityProblem_new_with_data(3, NC, NM_create_with_data(NM_DENSE, N, N, M), q,
                                              model);
  } else {
    problem = plasticityProblem_new_with_data(3, NC, NM_create_with_data(NM_DENSE, N, N, M), q,
                                              NULL);
    problem->model_type = PLASTICITY_MODEL_VON_MISES;
    problem->model.von_mises = plasticity_VonMises_model_new(p1);
    free(p2);
    free(p3);
  }
  return problem;
}

/* r_e (with the cohesion shift), v_e (modified strain rate), ||F_e||^2 */
static void element(PlasticityProblem* problem, int e, const double* z, const double* w,
                    double r[3], double v[3], double* f2) {
  double tmp[3];
  memcpy(r, z + 3 * e, sizeof(double) * 3);
  memcpy(v, w + 3 * e, sizeof(double) * 3);
  *f2 = 0.;
  if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
    double eta = problem->model.drucker_prager->eta[e];
    double theta = problem->model.drucker_prager->theta[e];
    r[0] += problem->model.drucker_prager->cohesion[e] / eta;
    plasticity_2d_unitary_compute_and_add_error(r, v, eta, theta, f2, tmp);
    v[0] += theta * hypot(v[1], v[2]);
  } else {
    plasticity_2d_vonmises_unitary_compute_and_add_error(
        r, v, problem->model.von_mises->sigma_y[e], f2, tmp);
  }
}

static double safe(double s) { return s > 1. ? s : 1.; }

static double error_with_rho(PlasticityProblem* problem, const double* z, int norm, int scope,
                             int reduction, int rho) {
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
  options->iparam[PLASTICITY_SCALING_RHO] = rho;
  options->iparam[PLASTICITY_SCALING_NORM] = norm;
  options->iparam[PLASTICITY_SCALING_SCOPE] = scope;
  options->iparam[PLASTICITY_SCALING_REDUCTION] = reduction;
  int n = 3 * problem->numberOfCones;
  double* w = (double*)malloc(n * sizeof(double));
  double* zc = (double*)malloc(n * sizeof(double));
  memcpy(zc, z, n * sizeof(double));
  double error;
  plasticity_2d_compute_error(problem, zc, w, 0., options, 0., &error);
  solver_options_delete(options);
  free(w);
  free(zc);
  return error;
}

/* rho = 1 (natural map of fc3d) */
static double error_with(PlasticityProblem* problem, const double* z, int norm, int scope,
                         int reduction) {
  return error_with_rho(problem, z, norm, scope, reduction, PLASTICITY_SCALING_RHO_ONE);
}

static void test_formulas(int drucker_prager) {
  printf("=== scaled errors, %s ===\n", drucker_prager ? "Drucker-Prager" : "Von Mises");
  PlasticityProblem* problem = random_problem(drucker_prager, 1.0);
  double z[N], w[N];
  for (int i = 0; i < N; i++) z[i] = 4. * uniform();
  plasticity_2d_compute_error(problem, z, w, 0., NULL, 0., &(double){0.}); /* w = M z + q */
  double error_default =
      error_with(problem, z, PLASTICITY_SCALING_EUCLIDEAN, PLASTICITY_SCALING_SCOPE_PER_CONE,
                 PLASTICITY_SCALING_REDUCTION_MAX);
  double abs2 = 0., max2 = 0., mixed2 = 0., nr = 0., nv = 0.;
  for (int e = 0; e < NC; e++) {
    double r[3], v[3], f2;
    element(problem, e, z, w, r, v, &f2);
    double a = safe(norm3(r)), b = safe(norm3(v));
    abs2 += f2;
    max2 += f2 / (fmin(a, b) * fmin(a, b));
    mixed2 += f2 / (safe(norm3(r) + norm3(v)) * safe(norm3(r) + norm3(v)));
    nr += norm3(r) * norm3(r);
    nv += norm3(v) * norm3(v);
  }
  double e_abs = error_with(problem, z, PLASTICITY_SCALING_NO, 0, 0);
  check(fabs(e_abs - sqrt(abs2)) <= 1e-14 * sqrt(abs2), "no scaling: unitary errors",
        fabs(e_abs - sqrt(abs2)));
  check(fabs(error_default - sqrt(max2)) <= 1e-14 * sqrt(max2),
        "default norm/scope/reduction: per cone, min of the safeguarded norms",
        fabs(error_default - sqrt(max2)));
  double e_mixed =
      error_with(problem, z, PLASTICITY_SCALING_EUCLIDEAN, PLASTICITY_SCALING_SCOPE_PER_CONE,
                 PLASTICITY_SCALING_REDUCTION_MIXED);
  check(fabs(e_mixed - sqrt(mixed2)) <= 1e-14 * sqrt(mixed2), "per cone, MIXED",
        fabs(e_mixed - sqrt(mixed2)));
  double e_vec =
      error_with(problem, z, PLASTICITY_SCALING_EUCLIDEAN, PLASTICITY_SCALING_SCOPE_PER_VECTOR,
                 PLASTICITY_SCALING_REDUCTION_MIXED);
  double expected = sqrt(abs2) / safe(sqrt(nr) + sqrt(nv));
  check(fabs(e_vec - expected) <= 1e-14 * expected, "per vector, MIXED",
        fabs(e_vec - expected));
  plasticityProblem_free(problem);
  printf("\n");
}

static void test_homogeneity(int drucker_prager) {
  printf("=== homogeneity (M, lambda q), lambda = 1e4, %s ===\n",
         drucker_prager ? "Drucker-Prager" : "Von Mises");
  const double lambda = 1e4;
  PlasticityProblem* p1 = random_problem(drucker_prager, 1.0);
  PlasticityProblem* pl = random_problem(drucker_prager, lambda);
  double z[N], zl[N];
  for (int i = 0; i < N; i++) {
    z[i] = 4. * uniform() + 5.; /* norms > 1 */
    zl[i] = lambda * z[i];
  }
  double a1 = error_with(p1, z, PLASTICITY_SCALING_NO, 0, 0);
  double al = error_with(pl, zl, PLASTICITY_SCALING_NO, 0, 0);
  check(fabs(al / a1 - lambda) <= 1e-10 * lambda, "absolute error multiplied by lambda",
        al / a1);
  double m1 = error_with(p1, z, 0, 0, PLASTICITY_SCALING_REDUCTION_MIXED);
  double ml = error_with(pl, zl, 0, 0, PLASTICITY_SCALING_REDUCTION_MIXED);
  check(fabs(ml / m1 - 1.) <= 1e-10, "MIXED scaled error unchanged", ml / m1);
  plasticityProblem_free(p1);
  plasticityProblem_free(pl);
  printf("\n");
}

/* rho_e = 1/||M_ee||_2: with (kappa M, kappa q) (strain rate in other units, same stress),
 * the error is unchanged; with rho = 1 it is not. With (M, lambda q) and lambda sigma, the
 * default error (MAX, rho_e) is unchanged when the norms are > 1. */
static void test_rho(int drucker_prager) {
  printf("=== rho_e = 1/||M_ee||: units of the strain rate, %s ===\n",
         drucker_prager ? "Drucker-Prager" : "Von Mises");
  const double kappa = 1e-6, lambda = 1e4;
  PlasticityProblem* p1 = random_problem(drucker_prager, 1.0);
  PlasticityProblem* pk = random_problem(drucker_prager, 1.0);
  PlasticityProblem* pl = random_problem(drucker_prager, lambda);
  NM_scal(kappa, pk->M);
  for (int i = 0; i < N; i++) pk->q[i] *= kappa;
  double z[N], zl[N];
  for (int i = 0; i < N; i++) {
    z[i] = 4. * uniform() + 5.;
    zl[i] = lambda * z[i];
  }
  double e1 = error_with_rho(p1, z, 0, 0, 0, PLASTICITY_SCALING_RHO_ELEMENT);
  double ek = error_with_rho(pk, z, 0, 0, 0, PLASTICITY_SCALING_RHO_ELEMENT);
  check(fabs(ek / e1 - 1.) < 1e-10, "(kappa M, kappa q): error unchanged (rho_e)", ek / e1);
  double o1 = error_with_rho(p1, z, 0, 0, 0, PLASTICITY_SCALING_RHO_ONE);
  double ok = error_with_rho(pk, z, 0, 0, 0, PLASTICITY_SCALING_RHO_ONE);
  check(fabs(ok / o1 - 1.) > 0.1, "(kappa M, kappa q): error changed with rho = 1", ok / o1);
  double g1 = error_with_rho(p1, z, 0, 0, 0, PLASTICITY_SCALING_RHO_GLOBAL);
  double gk = error_with_rho(pk, z, 0, 0, 0, PLASTICITY_SCALING_RHO_GLOBAL);
  check(fabs(gk / g1 - 1.) < 1e-10, "(kappa M, kappa q): error unchanged (global rho)",
        gk / g1);
  double el = error_with_rho(pl, zl, 0, 0, 0, PLASTICITY_SCALING_RHO_ELEMENT);
  check(fabs(el / e1 - 1.) < 1e-10,
        "(M, lambda q), lambda sigma: error unchanged (rho_e, MAX)", el / e1);
  /* options NULL = defaults of the NSGS (per cone, MIXED, rho_e) */
  double w[N], zc[N], enull;
  memcpy(zc, z, sizeof zc);
  plasticity_2d_compute_error(p1, zc, w, 0., NULL, 0., &enull);
  double emixed =
      error_with_rho(p1, z, PLASTICITY_SCALING_EUCLIDEAN, PLASTICITY_SCALING_SCOPE_PER_CONE,
                     PLASTICITY_SCALING_REDUCTION_MIXED, PLASTICITY_SCALING_RHO_ELEMENT);
  check(enull == emixed, "options NULL: per cone, MIXED, rho_e", enull - emixed);
  SolverOptions* nsgs = solver_options_create(PLASTICITY_2D_NSGS);
  double ensgs;
  memcpy(zc, z, sizeof zc);
  plasticity_2d_compute_error(p1, zc, w, 0., nsgs, 0., &ensgs);
  check(ensgs == emixed, "NSGS default options: per cone, MIXED, rho_e", ensgs - emixed);
  solver_options_delete(nsgs);
  plasticityProblem_free(p1);
  plasticityProblem_free(pk);
  plasticityProblem_free(pl);
  printf("\n");
}

static void test_footing(void) {
  printf("=== footing problem, rho_e and MIXED scaling, tolerance 1e-13 ===\n");
  PlasticityProblem* problem =
      plasticity_new_from_filename("./data/plasticity_2d_footing_100.dat");
  if (!problem) {
    check(0, "cannot read the data file", 0.);
    return;
  }
  int n = 3 * problem->numberOfCones;
  double* s = (double*)calloc(n, sizeof(double));
  double* w = (double*)calloc(n, sizeof(double));
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
  solver_options_update_internal(options, 0, PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-13;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
  options->iparam[PLASTICITY_SCALING_REDUCTION] = PLASTICITY_SCALING_REDUCTION_MIXED;
  options->iparam[PLASTICITY_SCALING_RHO] = PLASTICITY_SCALING_RHO_ELEMENT;
  int info = plasticity_2d_driver(problem, s, w, options);
  printf("  NSGS: info %d, %d iterations, scaled error %.2e\n", info,
         options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  check(info == 0, "converged", options->dparam[SICONOS_DPARAM_RESIDU]);
  double absolute = error_with(problem, s, PLASTICITY_SCALING_NO, 0, 0);
  printf("  absolute error %.2e\n", absolute);
  check(absolute < 1e-11, "absolute error at round-off", absolute);
  solver_options_delete(options);
  free(s);
  free(w);
  plasticityProblem_free(problem);
  printf("\n");
}

int main(void) {
  test_formulas(1);
  test_formulas(0);
  test_homogeneity(1);
  test_homogeneity(0);
  test_rho(1);
  test_rho(0);
  test_footing();
  printf("%d check(s) failed\n", failures);
  return failures ? 1 : 0;
}
