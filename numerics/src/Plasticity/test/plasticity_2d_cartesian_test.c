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
 * Cartesian stress coordinates (PlasticityProblem::stress_coordinates).
 *
 * The map between Cartesian and cone coordinates is written again here from the
 * definitions (sigma_N = -(sxx + syy)/sqrt2, sigma_T1 = (sxx - syy)/sqrt2,
 * sigma_T2 = sqrt2 sxy; power sxx exx + syy eyy + 2 sxy exy), independently of
 * plasticity_2d_stress_coordinates.c.
 *
 *  1. analytic: Von Mises with M = I, pure shear (sxy = sigma_y/sqrt3) and plane strain
 *     uniaxial loading (sxx - syy = 2 sigma_y/sqrt3), Voigt and Mandel;
 *  2. random coupled problems (Drucker-Prager with dilatancy and cohesion, Von Mises), given
 *     in cone coordinates and rewritten in Cartesian coordinates (Mandel, Voigt; dense,
 *     sparse-block and sparse storage): the Cartesian solution is the cone solution mapped
 *     back, the power is preserved and plasticity_2d_compute_error accepts it;
 *  3. footing problem of the Plasticity tests rewritten in Voigt coordinates.
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
#include "SparseBlockMatrix.h"
#include "plasticity_2d_compute_error.h"
#include "plasticity_2d_solvers.h"
#include "plasticity_2d_stress_coordinates.h"

static int failures = 0;

static void check(int ok, const char* what, double value) {
  printf("  %-64s %12.4e  %s\n", what, value, ok ? "ok" : "FAILURE");
  if (!ok) failures++;
}

static double dist(int n, const double* x, const double* y) {
  double s = 0.;
  for (int i = 0; i < n; i++) s += (x[i] - y[i]) * (x[i] - y[i]);
  return sqrt(s);
}

static const char* name(PlasticityStressCoordinates c) {
  return c == PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT ? "Voigt" : "Mandel";
}

/* shear component of the Cartesian stress / strain rate vectors */
static double shear_stress_factor(PlasticityStressCoordinates c) {
  return c == PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT ? 1.0 : M_SQRT2;
}
static double shear_strain_factor(PlasticityStressCoordinates c) {
  return c == PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT ? 2.0 : M_SQRT2;
}

/* P (row-major) such that sigma_c = P s, from the definitions */
static void matrix_P(PlasticityStressCoordinates c, double P[9]) {
  double f = shear_stress_factor(c); /* s[2] = f sxy */
  const double h = M_SQRT1_2;
  double A[9] = {-h, -h, 0., h, -h, 0., 0., 0., M_SQRT2 / f};
  memcpy(P, A, sizeof A);
}

/* s = P^{-1} sigma_c */
static void cone_to_cartesian_stress(PlasticityStressCoordinates c, int nc, const double* sc,
                                     double* s) {
  double f = shear_stress_factor(c);
  for (int e = 0; e < nc; e++) {
    double sN = sc[3 * e], s1 = sc[3 * e + 1], s2 = sc[3 * e + 2];
    double p = -sN / M_SQRT2, d = s1 / M_SQRT2; /* (sxx + syy)/2, (sxx - syy)/2 */
    s[3 * e] = p + d;
    s[3 * e + 1] = p - d;
    s[3 * e + 2] = f * s2 / M_SQRT2;
  }
}

/* w_s = P^T w_c */
static void cone_to_cartesian_strainrate(PlasticityStressCoordinates c, int nc,
                                         const double* wc, double* ws) {
  double P[9];
  matrix_P(c, P);
  for (int e = 0; e < nc; e++)
    for (int i = 0; i < 3; i++) {
      ws[3 * e + i] = 0.;
      for (int k = 0; k < 3; k++) ws[3 * e + i] += P[3 * k + i] * wc[3 * e + k];
    }
}

/* Cartesian problem: M_s = P^T M_c P (dense, column-major), q_s = P^T q_c */
static double* cartesian_matrix(PlasticityStressCoordinates c, int n, const double* Mc) {
  double P[9];
  matrix_P(c, P);
  double* Ms = (double*)calloc((size_t)(n * n), sizeof(double));
  for (int r = 0; r < n; r += 3)
    for (int c0 = 0; c0 < n; c0 += 3)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          double v = 0.;
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
              v += P[3 * k + i] * Mc[(r + k) + (c0 + l) * n] * P[3 * l + j];
          Ms[(r + i) + (c0 + j) * n] = v;
        }
  return Ms;
}

static NumericsMatrix* with_storage(int storage, int n, double* dense) {
  NumericsMatrix* D = NM_create_with_data(NM_DENSE, n, n, dense);
  if (storage == NM_DENSE) return D;
  NumericsMatrix* M;
  if (storage == NM_SPARSE)
    M = NM_copy_to_sparse(D, 0.);
  else {
    M = NM_create(NM_SPARSE_BLOCK, n, n);
    SBM_from_dense(3, (size_t)n, (size_t)n, dense, M->sparse_block);
  }
  NM_free(D);
  return M;
}

static int solve(PlasticityProblem* problem, double* sigma, double* w, double tol, int itermax,
                 int verbose_result) {
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = tol;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = itermax;
  int info = plasticity_2d_driver(problem, sigma, w, options);
  if (verbose_result)
    printf("  NSGS: info = %d, iterations = %d, error = %.3e\n", info,
           options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  solver_options_delete(options);
  return info;
}

/* ------------------------------------------------------------------------------------------
 */
static void test_von_mises_analytic(PlasticityStressCoordinates c) {
  printf("=== Von Mises, M = I, %s: pure shear and plane strain uniaxial loading ===\n",
         name(c));
  double sigma_y = 3.0, tau = 10.0;
  double* M = (double*)calloc(36, sizeof(double));
  for (int i = 0; i < 6; i++) M[i + 6 * i] = 1.;
  double* q = (double*)calloc(6, sizeof(double));
  q[2] = -tau; /* element 0: shear loading */
  q[3] = -tau; /* element 1: sxx loading */
  double* sy = (double*)malloc(2 * sizeof(double));
  sy[0] = sy[1] = sigma_y;
  PlasticityProblem* problem =
      plasticityProblem_new_with_data(3, 2, NM_create_with_data(NM_DENSE, 6, 6, M), q, NULL);
  problem->model_type = PLASTICITY_MODEL_VON_MISES;
  problem->model.von_mises = plasticity_VonMises_model_new(sy);
  problem->stress_coordinates = c;
  double s[6] = {0.}, w[6] = {0.};
  int info = solve(problem, s, w, 1e-14, 1000, 0);
  check(info == 0, "converged", 0.);
  double sxy = s[2] / shear_stress_factor(c);
  check(fabs(s[0]) + fabs(s[1]) < 1e-13 && fabs(sxy - sigma_y / sqrt(3.)) < 1e-13,
        "pure shear: sxy = sigma_y/sqrt(3)", fabs(sxy - sigma_y / sqrt(3.)));
  check(fabs(s[3] - s[4] - 2. * sigma_y / sqrt(3.)) < 1e-13,
        "plane strain uniaxial: sxx - syy = 2 sigma_y/sqrt(3)",
        fabs(s[3] - s[4] - 2. * sigma_y / sqrt(3.)));
  /* with the Frobenius metric (Mandel), the hydrostatic part is free: sxx + syy = tau */
  if (c == PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL)
    check(fabs(s[3] + s[4] - tau) < 1e-13, "plane strain uniaxial: sxx + syy = tau",
          fabs(s[3] + s[4] - tau));
  plasticityProblem_free(problem);
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 */
static double* random_cone_matrix(int nc) {
  int n = 3 * nc;
  double* M = (double*)calloc((size_t)(n * n), sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int same = (i / 3 == j / 3), next = abs(i / 3 - j / 3) == 1;
      if (same || next) M[i + j * n] = (same ? 0.3 : 0.1) * (2. * rand() / RAND_MAX - 1.);
    }
  for (int i = 0; i < n; i++) M[i + i * n] += 2.;
  return M;
}

static void test_random(int drucker_prager, PlasticityStressCoordinates c, int storage) {
  int nc = 4, n = 3 * nc;
  const char* storage_name =
      storage == NM_DENSE ? "dense" : (storage == NM_SPARSE ? "sparse" : "sparse block");
  printf("=== random %s problem, %s, %s ===\n",
         drucker_prager ? "Drucker-Prager" : "Von Mises", name(c), storage_name);
  srand(17 + drucker_prager);
  double* Mc = random_cone_matrix(nc);
  double* qc = (double*)malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) qc[i] = 2. * rand() / RAND_MAX - 1.;
  double* Ms = cartesian_matrix(c, n, Mc);
  double* qs = (double*)malloc((size_t)n * sizeof(double));
  cone_to_cartesian_strainrate(c, nc, qc, qs);

  PlasticityProblem* problems[2];
  for (int k = 0; k < 2; k++) {
    double* p1 = (double*)malloc((size_t)nc * sizeof(double));
    double* p2 = (double*)malloc((size_t)nc * sizeof(double));
    double* p3 = (double*)malloc((size_t)nc * sizeof(double));
    for (int e = 0; e < nc; e++) {
      p1[e] = drucker_prager ? 0.4 + 0.1 * e : 0.3 + 0.1 * e;
      p2[e] = 0.1 * (e % 2);
      p3[e] = 0.05 * e;
    }
    double* M = k ? Ms : Mc;
    double* q = k ? qs : qc;
    if (!k) {
      M = (double*)malloc((size_t)(n * n) * sizeof(double));
      memcpy(M, Mc, (size_t)(n * n) * sizeof(double));
      q = (double*)malloc((size_t)n * sizeof(double));
      memcpy(q, qc, (size_t)n * sizeof(double));
    }
    NumericsMatrix* NM =
        k ? with_storage(storage, n, M) : NM_create_with_data(NM_DENSE, n, n, M);
    if (drucker_prager) {
      Plasticity_DruckerPrager_model* model = plasticity_DruckerPrager_model_new(p1, p2);
      model->cohesion = p3;
      problems[k] = plasticityProblem_new_with_data(3, nc, NM, q, model);
    } else {
      problems[k] = plasticityProblem_new_with_data(3, nc, NM, q, NULL);
      problems[k]->model_type = PLASTICITY_MODEL_VON_MISES;
      problems[k]->model.von_mises = plasticity_VonMises_model_new(p1);
      free(p2);
      free(p3);
    }
  }
  problems[1]->stress_coordinates = c;

  double sc[12] = {0.}, wc[12] = {0.}, s[12] = {0.}, w[12] = {0.};
  int info0 = solve(problems[0], sc, wc, 1e-14, 10000, 0);
  int info1 = solve(problems[1], s, w, 1e-14, 10000, 1);
  check(info0 == 0 && info1 == 0, "both converged", 0.);

  /* number of elements on the yield surface (the test must not be purely elastic) */
  int plastic = 0;
  for (int e = 0; e < nc; e++) {
    double nT = hypot(sc[3 * e + 1], sc[3 * e + 2]);
    double radius = drucker_prager ? problems[0]->model.drucker_prager->eta[e] * sc[3 * e] +
                                         problems[0]->model.drucker_prager->cohesion[e]
                                   : sqrt(2. / 3.) * problems[0]->model.von_mises->sigma_y[e];
    if (fabs(nT - radius) < 1e-10) plastic++;
  }
  printf("  %d/%d elements on the yield surface\n", plastic, nc);
  check(plastic > 0, "some elements plastic", (double)plastic);

  double s_expected[12], w_expected[12];
  cone_to_cartesian_stress(c, nc, sc, s_expected);
  cone_to_cartesian_strainrate(c, nc, wc, w_expected);
  check(dist(n, s, s_expected) < 1e-12 * (1. + cblas_dnrm2(n, s, 1)),
        "Cartesian stress = cone solution mapped back", dist(n, s, s_expected));
  check(dist(n, w, w_expected) < 1e-12 * (1. + cblas_dnrm2(n, w, 1)),
        "Cartesian strain rate = cone solution mapped back", dist(n, w, w_expected));
  double power_c = cblas_ddot(n, sc, 1, wc, 1), power_s = cblas_ddot(n, s, 1, w, 1);
  check(fabs(power_c - power_s) < 1e-13 * (1. + fabs(power_c)), "power s.w_s = sigma.w",
        fabs(power_c - power_s));
  double error;
  plasticity_2d_compute_error(problems[1], s, w, 1e-12, NULL, 0., &error);
  check(error < 1e-12, "plasticity_2d_compute_error in Cartesian coordinates", error);

  /* the cone problem built by the library is the reference cone problem */
  PlasticityProblem* cone = plasticity_2d_problem_to_cone(problems[1]);
  NumericsMatrix* Mcone = NM_create(NM_DENSE, n, n);
  NM_to_dense(cone->M, Mcone);
  check(dist(n * n, Mcone->dense, Mc) < 1e-14 && dist(n, cone->q, qc) < 1e-14 &&
            cone->stress_coordinates == PLASTICITY_STRESS_COORDINATES_CONE,
        "plasticity_2d_problem_to_cone", dist(n * n, Mcone->dense, Mc));
  NM_free(Mcone);
  plasticityProblem_free(cone);

  plasticityProblem_free(problems[0]);
  plasticityProblem_free(problems[1]);
  free(Mc);
  free(qc);
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 */
static void test_footing(void) {
  printf("=== footing problem (100 cones) in Voigt coordinates ===\n");
  PlasticityProblem* problem =
      plasticity_new_from_filename("./data/plasticity_2d_footing_100.dat");
  if (!problem) {
    check(0, "cannot read the data file", 0.);
    return;
  }
  PlasticityStressCoordinates c = PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT;
  int nc = problem->numberOfCones, n = 3 * nc;
  PlasticityProblem* cartesian = plasticity_copy(problem);
  double* Ms = cartesian_matrix(c, n, problem->M->dense);
  NM_free(cartesian->M);
  cartesian->M = NM_create_with_data(NM_DENSE, n, n, Ms);
  cone_to_cartesian_strainrate(c, nc, problem->q, cartesian->q);
  cartesian->stress_coordinates = c;

  double* sc = (double*)calloc(n, sizeof(double));
  double* wc = (double*)calloc(n, sizeof(double));
  double* s = (double*)calloc(n, sizeof(double));
  double* w = (double*)calloc(n, sizeof(double));
  double* s_expected = (double*)calloc(n, sizeof(double));
  int info0 = solve(problem, sc, wc, 1e-12, 100000, 1);
  int info1 = solve(cartesian, s, w, 1e-12, 100000, 1);
  check(info0 == 0 && info1 == 0, "both converged", 0.);
  cone_to_cartesian_stress(c, nc, sc, s_expected);
  double d = dist(n, s, s_expected) / (1. + cblas_dnrm2(n, s_expected, 1));
  check(d < 1e-8, "Cartesian solution = cone solution mapped back (relative)", d);
  double error;
  plasticity_2d_compute_error(cartesian, s, w, 1e-10, NULL, 0., &error);
  check(error < 1e-10, "plasticity_2d_compute_error in Cartesian coordinates", error);

  free(sc);
  free(wc);
  free(s);
  free(w);
  free(s_expected);
  plasticityProblem_free(cartesian);
  plasticityProblem_free(problem);
  printf("\n");
}

int main(void) {
  test_von_mises_analytic(PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT);
  test_von_mises_analytic(PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL);
  int storages[3] = {NM_DENSE, NM_SPARSE_BLOCK, NM_SPARSE};
  for (int dp = 0; dp < 2; dp++)
    for (int k = 0; k < 3; k++) {
      test_random(dp, PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL, storages[k]);
      test_random(dp, PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT, storages[k]);
    }
  test_footing();
  printf("%d check(s) failed\n", failures);
  return failures ? 1 : 0;
}
