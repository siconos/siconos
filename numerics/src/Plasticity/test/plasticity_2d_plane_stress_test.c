// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*
 * Plane stress Von Mises (PLASTICITY_VONMISES_PLANE_STRESS): yield set
 *   f(s) = sxx^2 + syy^2 - sxx syy + 3 sxy^2 <= sigma_y^2,
 * rewritten by plasticity_2d_driver as a ball problem (PLASTICITY_VONMISES_BALL).
 *
 *  1. ball local solvers: return mapping, fixed point and semi-smooth Newton on 1000
 *     pseudo-random non-symmetric local problems;
 *  2. manufactured solutions: for target stresses s (uniaxial tension and compression,
 *     equibiaxial, pure shear, random points of the yield surface, an interior point) and a
 *     coupled non-symmetric positive definite M, q = -M s - lambda grad f(s) (associated flow
 *     rule, written here in Cartesian Voigt / Mandel and in cone coordinates, independently
 *     of the library). The solution is unique and must be s, with all the local solvers and
 *     storages;
 *  3. plasticity_2d_problem_change_of_variables: the ball problem gives the same solution.
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
#include "plasticity_2d_vonmises.h"

static int failures = 0;

static void check(int ok, const char* what, double value) {
  printf("  %-66s %12.4e  %s\n", what, value, ok ? "ok" : "FAILURE");
  if (!ok) failures++;
}

static double dist(int n, const double* x, const double* y) {
  double s = 0.;
  for (int i = 0; i < n; i++) s += (x[i] - y[i]) * (x[i] - y[i]);
  return sqrt(s);
}

static double uniform(void) { return 2. * rand() / RAND_MAX - 1.; }

/* ------------------------------------------------------------------------------------------
 * 1. ball local solvers
 */
static double ball_residual(const double M[9], const double q[3], double r,
                            const double s[3]) {
  double y[3];
  for (int i = 0; i < 3; i++)
    y[i] = s[i] - (q[i] + M[i] * s[0] + M[i + 3] * s[1] + M[i + 6] * s[2]);
  plasticity_2d_projectionOnVonMisesBall(y, r);
  return dist(3, s, y);
}

static void test_ball_local_solvers(void) {
  printf("=== ball local solvers, 1000 random non-symmetric local problems ===\n");
  srand(3);
  double max_rm = 0., max_fp = 0., max_nsn = 0.;
  int plastic = 0, fails = 0;
  long it_fp = 0, it_nsn = 0;
  for (int k = 0; k < 1000; k++) {
    double M[9], q[3], r = 0.5 + 0.5 * (rand() / (double)RAND_MAX);
    for (int i = 0; i < 9; i++) M[i] = 0.3 * uniform();
    M[0] += 2.;
    M[4] += 1. + rand() / (double)RAND_MAX;
    M[8] += 1.5;
    for (int i = 0; i < 3; i++) q[i] = 3. * uniform();
    double s_rm[3], s_fp[3] = {0.}, s_nsn[3] = {0.};
    fails += plasticity_2d_vonMises_ball_return_mapping(M, q, r, s_rm) == 1;
    double n = sqrt(s_rm[0] * s_rm[0] + s_rm[1] * s_rm[1] + s_rm[2] * s_rm[2]);
    if (fabs(n - r) < 1e-12) plastic++;
    max_rm = fmax(max_rm, ball_residual(M, q, r, s_rm));

    Plasticity2DVonMisesLocalParameters p = {
        .max_iter = 2000,
        .tol = 1e-14,
        .rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE,
        .linesearch_max_iter = 10,
        .ball = 1};
    fails += plasticity_2d_vonMises_fixed_point(M, q, r, s_fp, &p);
    it_fp += p.iter_done;
    max_fp = fmax(max_fp, dist(3, s_fp, s_rm));

    Plasticity2DVonMisesLocalParameters pn = {
        .max_iter = 50,
        .tol = 1e-14,
        .rho_strategy = PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM,
        .linesearch = PLASTICITY_NSN_LINESEARCH_ARMIJO,
        .linesearch_max_iter = 20,
        .hybrid = 0,
        .ball = 1};
    fails += plasticity_2d_vonMises_semismooth_newton(M, q, r, s_nsn, &pn);
    it_nsn += pn.iter_done;
    max_nsn = fmax(max_nsn, dist(3, s_nsn, s_rm));
  }
  printf("  %d/1000 plastic, mean iterations: fixed point %.1f, Newton %.1f\n", plastic,
         it_fp / 1000., it_nsn / 1000.);
  check(fails == 0, "all solvers succeeded", (double)fails);
  check(plastic > 300 && plastic < 1000, "elastic and plastic cases", (double)plastic);
  check(max_rm < 1e-14, "return mapping: natural-map residual", max_rm);
  check(max_fp < 1e-12, "fixed point = return mapping", max_fp);
  check(max_nsn < 1e-12, "semi-smooth Newton = return mapping", max_nsn);
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 * 2. manufactured solutions
 */
#define NC 8

/* targets (sxx, syy, sxy) for sigma_y = 1 (scaled per element), lambda = 0: interior */
static void targets(double s[NC][3], double lambda[NC]) {
  const double t = 1. / sqrt(3.);
  double base[6][3] = {{1., 0., 0.},  /* uniaxial tension */
                       {0., -1., 0.}, /* uniaxial compression */
                       {1., 1., 0.},  /* equibiaxial tension */
                       {0., 0., t},   /* pure shear */
                       {0., 0., 0.},  /* random points of the yield surface: below */
                       {0., 0., 0.}};
  for (int e = 0; e < 6; e++) memcpy(s[e], base[e], sizeof base[e]);
  for (int e = 4; e < 6; e++) {
    double a = uniform(), b = uniform(), c = uniform();
    double f = sqrt(a * a + b * b - a * b + 3. * c * c);
    s[e][0] = a / f;
    s[e][1] = b / f;
    s[e][2] = c / f;
  }
  /* interior point */
  s[6][0] = 0.3;
  s[6][1] = -0.2;
  s[6][2] = 0.1;
  /* equibiaxial compression */
  s[7][0] = s[7][1] = -1.;
  s[7][2] = 0.;
  for (int e = 0; e < NC; e++) lambda[e] = (e == 6) ? 0. : 0.4 + 0.1 * e;
}

static double plane_stress_f(const double c[3]) { /* c = (sxx, syy, sxy) */
  return sqrt(c[0] * c[0] + c[1] * c[1] - c[0] * c[1] + 3. * c[2] * c[2]);
}

/* stress vector and gradient of f^2/2 in the given coordinates, from (sxx, syy, sxy) */
static void stress_and_gradient(PlasticityStressCoordinates coords, const double c[3],
                                double s[3], double g[3]) {
  double sxx = c[0], syy = c[1], sxy = c[2];
  switch (coords) {
    case PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT: /* f^2 = ... + 3 s3^2 */
      s[0] = sxx, s[1] = syy, s[2] = sxy;
      g[0] = sxx - 0.5 * syy, g[1] = syy - 0.5 * sxx, g[2] = 3. * sxy;
      break;
    case PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL: /* f^2 = ... + 3/2 s3^2 */
      s[0] = sxx, s[1] = syy, s[2] = M_SQRT2 * sxy;
      g[0] = sxx - 0.5 * syy, g[1] = syy - 0.5 * sxx, g[2] = 1.5 * s[2];
      break;
    default: /* cone: f^2 = sN^2/2 + 3/2 |sT|^2 */
      s[0] = -(sxx + syy) / M_SQRT2, s[1] = (sxx - syy) / M_SQRT2, s[2] = M_SQRT2 * sxy;
      g[0] = 0.5 * s[0], g[1] = 1.5 * s[1], g[2] = 1.5 * s[2];
  }
}

/* (sxx, syy, sxy) from the stress vector */
static void cartesian(PlasticityStressCoordinates coords, const double s[3], double c[3]) {
  switch (coords) {
    case PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT:
      c[0] = s[0], c[1] = s[1], c[2] = s[2];
      break;
    case PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL:
      c[0] = s[0], c[1] = s[1], c[2] = s[2] / M_SQRT2;
      break;
    default: {
      double p = -s[0] / M_SQRT2, d = s[1] / M_SQRT2;
      c[0] = p + d, c[1] = p - d, c[2] = s[2] / M_SQRT2;
    }
  }
}

static const char* coords_name(PlasticityStressCoordinates c) {
  return c == PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT
             ? "Voigt"
             : (c == PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL ? "Mandel" : "cone");
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

static PlasticityProblem* manufactured_problem(PlasticityStressCoordinates coords, int storage,
                                               double s_exact[3 * NC]) {
  int n = 3 * NC;
  srand(11);
  double c[NC][3], lambda[NC];
  targets(c, lambda);
  double* sigma_y = (double*)malloc(NC * sizeof(double));
  double* M = (double*)calloc((size_t)(n * n), sizeof(double));
  double* q = (double*)calloc((size_t)n, sizeof(double));
  /* coupled non-symmetric M with a positive definite symmetric part */
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      if (abs(i / 3 - j / 3) <= 1) M[i + j * n] = (i / 3 == j / 3 ? 0.3 : 0.1) * uniform();
  for (int i = 0; i < n; i++) M[i + i * n] += 2.;
  double g[3 * NC];
  for (int e = 0; e < NC; e++) {
    sigma_y[e] = 1. + 0.5 * e; /* scale the target */
    for (int i = 0; i < 3; i++) c[e][i] *= sigma_y[e];
    stress_and_gradient(coords, c[e], s_exact + 3 * e, g + 3 * e);
  }
  /* q = -M s - lambda grad (f^2/2), lambda >= 0 */
  for (int i = 0; i < n; i++) {
    q[i] = -lambda[i / 3] * g[i];
    for (int j = 0; j < n; j++) q[i] -= M[i + j * n] * s_exact[j];
  }
  PlasticityProblem* problem =
      plasticityProblem_new_with_data(3, NC, with_storage(storage, n, M), q, NULL);
  problem->model_type = PLASTICITY_MODEL_VON_MISES;
  problem->model.von_mises = plasticity_VonMises_model_new(sigma_y);
  problem->model.von_mises->hypothesis = PLASTICITY_VONMISES_PLANE_STRESS;
  problem->stress_coordinates = coords;
  return problem;
}

static void test_manufactured(PlasticityStressCoordinates coords, int storage, int local_id) {
  const char* storage_name =
      storage == NM_DENSE ? "dense" : (storage == NM_SPARSE ? "sparse" : "sparse block");
  printf("=== manufactured plane stress solutions, %s, %s, local solver %s ===\n",
         coords_name(coords), storage_name,
         local_id ? solver_options_id_to_name(local_id) : "default");
  int n = 3 * NC;
  double s_exact[3 * NC];
  PlasticityProblem* problem = manufactured_problem(coords, storage, s_exact);
  double s[3 * NC] = {0.}, w[3 * NC] = {0.};
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
  if (local_id) solver_options_update_internal(options, 0, local_id);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
  int info = plasticity_2d_driver(problem, s, w, options);
  printf("  NSGS: info = %d, iterations = %d, error = %.3e\n", info,
         options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  solver_options_delete(options);
  check(info == 0, "converged", 0.);
  double d = dist(n, s, s_exact) / cblas_dnrm2(n, s_exact, 1);
  check(d < 1e-11, "solution = manufactured stress (relative)", d);
  double max_f = 0.;
  for (int e = 0; e < NC; e++) {
    if (e == 6) continue; /* interior */
    double ce[3];
    cartesian(coords, s + 3 * e, ce);
    max_f = fmax(max_f, fabs(plane_stress_f(ce) / problem->model.von_mises->sigma_y[e] - 1.));
  }
  check(max_f < 1e-11, "plastic elements: sigma_eq = sigma_y (plane stress formula)", max_f);
  double error;
  plasticity_2d_compute_error(problem, s, w, 1e-12, NULL, 0., &error);
  check(error < 1e-11, "plasticity_2d_compute_error (plane stress)", error);
  check(problem->model.von_mises->hypothesis == PLASTICITY_VONMISES_PLANE_STRESS,
        "problem unchanged (hypothesis)", 0.);
  plasticityProblem_free(problem);
  printf("\n");
}

/* ------------------------------------------------------------------------------------------
 * 3. the ball problem built by plasticity_2d_problem_change_of_variables
 */
static void test_change_of_variables(void) {
  printf("=== plasticity_2d_problem_change_of_variables (Voigt, plane stress -> ball) ===\n");
  int n = 3 * NC;
  double s_exact[3 * NC];
  PlasticityProblem* problem =
      manufactured_problem(PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT, NM_DENSE, s_exact);
  PlasticityProblem* ball = plasticity_2d_problem_change_of_variables(problem);
  check(ball->model.von_mises->hypothesis == PLASTICITY_VONMISES_BALL &&
            ball->stress_coordinates == PLASTICITY_STRESS_COORDINATES_CONE,
        "ball problem in cone coordinates", 0.);
  double y[3 * NC] = {0.}, w[3 * NC] = {0.}, s[3 * NC];
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  int info = plasticity_2d_driver(ball, y, w, options);
  solver_options_delete(options);
  double T[9], Tinv[9];
  plasticity_2d_change_of_variables(problem, T, Tinv);
  plasticity_2d_block_map(T, 0, NC, y, s);
  double d = dist(n, s, s_exact) / cblas_dnrm2(n, s_exact, 1);
  check(info == 0 && d < 1e-11, "s = T y", d);
  int on_ball = 0;
  for (int e = 0; e < NC; e++)
    on_ball += fabs(sqrt(y[3 * e] * y[3 * e] + y[3 * e + 1] * y[3 * e + 1] +
                         y[3 * e + 2] * y[3 * e + 2]) -
                    ball->model.von_mises->sigma_y[e]) < 1e-11;
  check(on_ball == NC - 1, "plastic elements on the ball ||y|| = sigma_y", (double)on_ball);
  plasticityProblem_free(ball);
  plasticityProblem_free(problem);
  printf("\n");
}

int main(void) {
  test_ball_local_solvers();
  PlasticityStressCoordinates coords[3] = {PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT,
                                           PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL,
                                           PLASTICITY_STRESS_COORDINATES_CONE};
  for (int k = 0; k < 3; k++) test_manufactured(coords[k], NM_DENSE, 0);
  test_manufactured(PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT, NM_SPARSE_BLOCK, 0);
  test_manufactured(PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT, NM_SPARSE, 0);
  int local_ids[3] = {PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING,
                      PLASTICITY_2D_ONECONE_VONMISES_FP,
                      PLASTICITY_2D_ONECONE_ProjectionOnCone};
  for (int k = 0; k < 3; k++)
    test_manufactured(PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT, NM_DENSE, local_ids[k]);
  test_change_of_variables();
  printf("%d check(s) failed\n", failures);
  return failures ? 1 : 0;
}
