// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*
 * Tests of the Mohr-Coulomb (MIXED_MODEL_DRUCKER_PRAGER) blocks of MixedProblem.
 *
 * Problem for one cone (stress sigma in R^3, w = M sigma + q):
 *
 *   K_eta = { sigma : ||sigma_T|| <= eta sigma_N },
 *   K_eta  ∋ sigma  ⊥  w + theta ||w_T|| e_N  ∈ K_eta^*,
 *
 * theta being the dilatancy coefficient (same convention as Plasticity/plasticity_2d_*).
 *
 * The solutions are checked independently of the MixedProblems code:
 *  - natural-map residual from plasticity_2d_unitary_compute_and_add_error (and
 *    fc2d_unitary_compute_and_add_error for friction blocks), with w recomputed as M z + q;
 *  - the w returned by the solver against M z + q;
 *  - analytic solutions (M = I) or the solution of the Plasticity NSGS solver.
 *
 * Matrices are dense, column-major (Siconos convention). Some test matrices are
 * non-symmetric on purpose, to detect transposed products.
 */

#include "NonSmoothDrivers.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "mixed_test_tools.h"

#ifndef PLASTICITY_DATA_DIR
#define PLASTICITY_DATA_DIR "."
#endif

/* Reference solution with the Plasticity NSGS solver (all cones Mohr-Coulomb). */
static int solve_plasticity_reference(int nc, const double* M, const double* q,
                                      const double* eta, const double* theta,
                                      const double* cohesion, double* z) {
  int n = 3 * nc;
  double* Mc = (double*)malloc((size_t)n * n * sizeof(double));
  double* qc = (double*)malloc(n * sizeof(double));
  double* etac = (double*)malloc(nc * sizeof(double));
  double* thetac = (double*)malloc(nc * sizeof(double));
  memcpy(Mc, M, (size_t)n * n * sizeof(double));
  memcpy(qc, q, n * sizeof(double));
  memcpy(etac, eta, nc * sizeof(double));
  memcpy(thetac, theta, nc * sizeof(double));
  Plasticity_DruckerPrager_model* model = plasticity_DruckerPrager_model_new(etac, thetac);
  if (cohesion) {
    model->cohesion = (double*)malloc(nc * sizeof(double));
    memcpy(model->cohesion, cohesion, nc * sizeof(double));
  }
  PlasticityProblem* P = plasticityProblem_new_with_data(
      3, nc, NM_create_with_data(NM_DENSE, n, n, Mc), qc, model);
  double* w = (double*)calloc(n, sizeof(double));
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
  int info = plasticity_2d_driver(P, z, w, options);
  printf("  reference (Plasticity NSGS): info = %d, iterations = %d, error = %.3e\n", info,
         options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  solver_options_delete(options);
  plasticityProblem_free(P);
  free(w);
  return info;
}

/* ------------------------------------------------------------------------------------------
 */
/* Tests */
/* ------------------------------------------------------------------------------------------
 */

/* M = I, theta = 0: sigma = P_K(-q) (Moreau decomposition). One cone per regime:
   stick (interior), separation (sigma = 0), sliding (boundary). */
static int test_identity_theta0(void) {
  printf("=== Mohr-Coulomb, M = I, theta = 0: stick / separation / sliding ===\n");
  int n = 9;
  double* M = (double*)calloc(n * n, sizeof(double));
  for (int i = 0; i < n; i++) M[i + i * n] = 1.;
  double q[9] = {-1., -0.1, 0.05, 1., 0.2, 0., -1., -1., 0.};
  double eta[3] = {0.5, 0.5, 0.5}, theta[3] = {0., 0., 0.};
  /* -q = (1, 1, 0) projected on K_0.5: n = (1 + 0.5)/(1 + 0.25) = 1.2, t = 0.5 n = 0.6 */
  double z_expected[9] = {1., 0.1, -0.05, 0., 0., 0., 1.2, 0.6, 0.};
  BlockSpec blocks[3];
  for (int b = 0; b < 3; b++)
    blocks[b] = (BlockSpec){MIXED_MODEL_DRUCKER_PRAGER, 3 * b, 3, 1, &eta[b], &theta[b]};

  double z[9] = {0.}, w[9] = {0.};
  int solver_info = solve_mixed(n, M, q, 3, blocks, z, w, NM_DENSE);
  int info = check_solution("identity_theta0", solver_info, n, M, q, 3, blocks, z, w,
                            z_expected, 1e-12);
  free(M);
  return info;
}

/* M = I, -q = (1, 1, 0), eta = theta = 0.5. Sliding with sigma = (n, eta n, 0):
   w = (n-1, eta n - 1, 0), w + theta|w_T| e_N orthogonal to sigma and on the boundary of K^*
   gives n = 1, i.e. sigma = (1, 0.5, 0), w = (0, -0.5, 0): no normal (dilatant) rate. */
static int test_identity_dilatancy(void) {
  printf("=== Mohr-Coulomb, M = I, eta = theta = 0.5, sliding ===\n");
  int n = 3;
  double M[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  double q[3] = {-1., -1., 0.};
  double eta[1] = {0.5}, theta[1] = {0.5};
  double z_expected[3] = {1., 0.5, 0.};
  BlockSpec blocks[1] = {{MIXED_MODEL_DRUCKER_PRAGER, 0, 3, 1, eta, theta}};
  double z[3] = {0.}, w[3] = {0.};
  int solver_info = solve_mixed(n, M, q, 1, blocks, z, w, NM_DENSE);
  return check_solution("identity_dilatancy", solver_info, n, M, q, 1, blocks, z, w,
                        z_expected, 1e-8);
}

/* Two coupled cones, non-symmetric M (diagonally dominant), with or without cohesion.
   Reference: Plasticity NSGS (plasticity_2d_driver, which translates the cohesion). */
static const double M2_rows[36] = {2.0, 0.3, 0.1, 0.2, 0.1, 0.0,  /* */
                                   0.1, 1.5, 0.2, 0.0, 0.2, 0.1,  /* */
                                   0.0, 0.1, 1.8, 0.1, 0.0, 0.2,  /* */
                                   0.3, 0.0, 0.1, 2.2, 0.2, 0.1,  /* */
                                   0.1, 0.2, 0.0, 0.0, 1.6, 0.3,  /* */
                                   0.0, 0.1, 0.2, 0.1, 0.1, 1.7}; /* */

static int test_coupled_nonsymmetric(int cones_per_block, int with_cohesion) {
  printf("=== Mohr-Coulomb, 2 coupled cones, non-symmetric M, %s%s ===\n",
         cones_per_block == 1 ? "one block per cone" : "one block with 2 cones",
         with_cohesion ? ", cohesion" : "");
  int n = 6, nc = 2;
  double* M = colmajor_from_rowmajor(n, M2_rows);
  double q[6] = {-1., 0.8, -0.3, -0.5, -0.6, 0.4};
  double eta[2] = {0.5, 0.4}, theta[2] = {0.3, 0.};
  double cohesion_values[2] = {0.2, 0.3};
  double* cohesion = with_cohesion ? cohesion_values : NULL;

  double z_ref[6] = {0.};
  int info = solve_plasticity_reference(nc, M, q, eta, theta, cohesion, z_ref);
  BlockSpec ref_block[1] = {{MIXED_MODEL_DRUCKER_PRAGER, 0, 3, nc, eta, theta, cohesion}};
  printf("  reference scaled residual = %.3e\n", mixed_residual(n, M, q, 1, ref_block, z_ref));

  BlockSpec blocks[2];
  int nb;
  if (cones_per_block == 1) {
    nb = 2;
    blocks[0] = (BlockSpec){MIXED_MODEL_DRUCKER_PRAGER,    0, 3, 1, &eta[0], &theta[0],
                            cohesion ? &cohesion[0] : NULL};
    blocks[1] = (BlockSpec){MIXED_MODEL_DRUCKER_PRAGER,    3, 3, 1, &eta[1], &theta[1],
                            cohesion ? &cohesion[1] : NULL};
  } else {
    nb = 1;
    blocks[0] = (BlockSpec){MIXED_MODEL_DRUCKER_PRAGER, 0, 3, 2, eta, theta, cohesion};
  }
  double z[6] = {0.}, w[6] = {0.};
  int solver_info = solve_mixed(n, M, q, nb, blocks, z, w, NM_DENSE);
  info += check_solution(cones_per_block == 1 ? "coupled_nonsymmetric_1cone_per_block"
                                              : "coupled_nonsymmetric_2cones_1block",
                         solver_info, n, M, q, nb, blocks, z, w, z_ref, 1e-8);
  free(M);
  return info;
}

/* Footing problem of the Plasticity tests (100 cones, dense SPD M ~ 1e-7, eta = theta = 0.1).
 * Solved with one block per cone, then with one block holding all the cones: both must have a
 * small (scaled) residual, agree, and agree with the Plasticity NSGS solution. */
static int test_footing(void) {
  printf("=== Mohr-Coulomb, footing problem (100 cones) ===\n");
  const char* filename = PLASTICITY_DATA_DIR "/plasticity_2d_footing_100.dat";
  PlasticityProblem* P = plasticity_new_from_filename(filename);
  if (!P) {
    printf("  cannot read %s\nTest footing: FAILED\n\n", filename);
    return 1;
  }
  int nc = P->numberOfCones, n = 3 * nc;
  double* Mdata = (double*)calloc((size_t)n * n, sizeof(double));
  NumericsMatrix* Md = NM_create_with_data(NM_DENSE, n, n, Mdata);
  NM_to_dense(P->M, Md);
  double* eta = P->model.drucker_prager->eta;
  double* theta = P->model.drucker_prager->theta;

  BlockSpec one_block[1] = {{MIXED_MODEL_DRUCKER_PRAGER, 0, 3, nc, eta, theta}};
  BlockSpec* per_cone = (BlockSpec*)malloc(nc * sizeof(BlockSpec));
  for (int c = 0; c < nc; c++)
    per_cone[c] = (BlockSpec){MIXED_MODEL_DRUCKER_PRAGER, 3 * c, 3, 1, &eta[c], &theta[c]};

  double* z1 = (double*)calloc(n, sizeof(double));
  double* z2 = (double*)calloc(n, sizeof(double));
  double* w = (double*)calloc(n, sizeof(double));

  printf("-- one block per cone\n");
  int solver_info = solve_mixed(n, Mdata, P->q, nc, per_cone, z1, w, NM_DENSE);
  int info = check_solution("footing_one_block_per_cone", solver_info, n, Mdata, P->q, nc,
                            per_cone, z1, w, NULL, 0.);

  printf("-- one block with all the cones (compared with the previous solution)\n");
  solver_info = solve_mixed(n, Mdata, P->q, 1, one_block, z2, w, NM_DENSE);
  info += check_solution("footing_one_block", solver_info, n, Mdata, P->q, 1, one_block, z2, w,
                         z1, 1e-8);

  printf("-- Plasticity NSGS solution\n");
  double* z_ref = (double*)calloc(n, sizeof(double));
  solve_plasticity_reference(nc, Mdata, P->q, eta, theta, NULL, z_ref);
  double d_ref = dist2(n, z_ref, z1) / (1. + norm2(n, z1));
  printf("  scaled residual = %.3e, |z_ref - z| / (1+|z|) = %.3e\n",
         mixed_residual(n, Mdata, P->q, 1, one_block, z_ref), d_ref);
  if (!(d_ref <= 1e-7)) { /* the mixed solve stops at 1e-10 (relative to |q|) */
    printf(
        "  FAILURE: mixed and Plasticity NSGS solutions differ\nTest footing_reference: "
        "FAILED\n\n");
    info++;
  } else
    printf("Test footing_reference: PASSED\n\n");

  free(z1);
  free(z2);
  free(w);
  free(z_ref);
  free(per_cone);
  NM_free(Md);
  plasticityProblem_free(P);
  return info;
}

/* Friction 2D contact coupled with a Mohr-Coulomb cone, non-symmetric M (5x5). */
static int test_friction2d_and_mohr_coulomb(void) {
  printf("=== Friction 2D + Mohr-Coulomb, coupled, non-symmetric M ===\n");
  int n = 5;
  const double M_rows[25] = {2.0, 0.2, 0.3, 0.1, 0.0,  /* */
                             0.1, 1.5, 0.0, 0.2, 0.1,  /* */
                             0.2, 0.1, 1.8, 0.3, 0.1,  /* */
                             0.0, 0.3, 0.1, 1.6, 0.2,  /* */
                             0.1, 0.0, 0.2, 0.1, 1.4}; /* */
  double* M = colmajor_from_rowmajor(n, M_rows);
  double q[5] = {-1., 0.9, -1., -0.8, 0.3};
  double mu[1] = {0.3}, eta[1] = {0.5}, theta[1] = {0.2};
  BlockSpec blocks[2] = {{MIXED_MODEL_FRICTION_2D, 0, 2, 1, mu, NULL},
                         {MIXED_MODEL_DRUCKER_PRAGER, 2, 3, 1, eta, theta}};
  double z[5] = {0.}, w[5] = {0.};
  int solver_info = solve_mixed(n, M, q, 2, blocks, z, w, NM_DENSE);
  int info = check_solution("friction2d_and_mohr_coulomb", solver_info, n, M, q, 2, blocks, z,
                            w, NULL, 0.);
  free(M);
  return info;
}

int main(void) {
  int failures = 0;
  failures += test_identity_theta0();
  failures += test_identity_dilatancy();
  failures += test_coupled_nonsymmetric(1, 0);
  failures += test_coupled_nonsymmetric(2, 0);
  failures += test_coupled_nonsymmetric(1, 1);
  failures += test_coupled_nonsymmetric(2, 1);
  failures += test_footing();
  failures += test_friction2d_and_mohr_coulomb();
  printf("%d test(s) FAILED\n", failures);
  return failures ? 1 : 0;
}
