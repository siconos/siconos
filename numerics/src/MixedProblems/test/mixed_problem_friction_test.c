// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*
 * Tests of the friction blocks of MixedProblem and of the storage of M:
 *  - friction 2D and 3D blocks with several contacts (one block vs one block per contact);
 *  - rolling friction 2D and 3D blocks;
 *  - blocks of different dimensions (generic extraction of the diagonal blocks);
 *  - dense and sparse (NM_SPARSE) storage of M.
 * The matrices are non-symmetric, diagonally dominant. Solutions are checked with the
 * natural-map residual of each element (unitary error functions of the single-model solvers),
 * w against M z + q, and against each other.
 */

#include "mixed_test_tools.h"

/* Non-symmetric diagonally dominant n x n matrix (dense, column-major). */
static double* test_matrix(int n) {
  double* M = (double*)calloc((size_t)n * n, sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      M[i + j * n] = (i == j) ? 2.0 + 0.1 * i : 0.15 * sin(7.0 * i + 3.0 * j + 1.0);
  return M;
}

/* Solve the same problem with two block decompositions (and/or storages), check both and
   compare the solutions. */
static int solve_and_compare(const char* name, int n, double* M, double* q, int nb1,
                             const BlockSpec* blocks1, NM_types storage1, int nb2,
                             const BlockSpec* blocks2, NM_types storage2) {
  char name1[256], name2[256];
  snprintf(name1, sizeof name1, "%s (1)", name);
  snprintf(name2, sizeof name2, "%s (2), compared with (1)", name);
  double* z1 = (double*)calloc(n, sizeof(double));
  double* z2 = (double*)calloc(n, sizeof(double));
  double* w = (double*)calloc(n, sizeof(double));
  int info = solve_mixed(n, M, q, nb1, blocks1, z1, w, storage1);
  int failures = check_solution(name1, info, n, M, q, nb1, blocks1, z1, w, NULL, 0.);
  info = solve_mixed(n, M, q, nb2, blocks2, z2, w, storage2);
  failures += check_solution(name2, info, n, M, q, nb2, blocks2, z2, w, z1, 1e-8);
  free(z1);
  free(z2);
  free(w);
  return failures;
}

/* Friction 2D, 3 contacts: one block with 3 contacts vs 3 blocks. */
static int test_friction2d_contacts(void) {
  printf("=== Friction 2D, 3 contacts: one block vs one block per contact ===\n");
  int n = 6;
  double* M = test_matrix(n);
  double q[6] = {-1., 0.5, -0.8, -1.2, 0.4, 0.1};
  double mu[3] = {0.3, 0.5, 0.2};
  BlockSpec one[1] = {{MIXED_MODEL_FRICTION_2D, 0, 2, 3, mu, NULL}};
  BlockSpec per_contact[3];
  for (int c = 0; c < 3; c++)
    per_contact[c] = (BlockSpec){MIXED_MODEL_FRICTION_2D, 2 * c, 2, 1, &mu[c], NULL};
  int failures = solve_and_compare("friction2d_contacts", n, M, q, 3, per_contact, NM_DENSE, 1,
                                   one, NM_DENSE);
  free(M);
  return failures;
}

/* Friction 3D, 3 contacts: one block vs 3 blocks, then dense vs sparse M. */
static int test_friction3d_contacts(void) {
  printf("=== Friction 3D, 3 contacts: one block vs one block per contact ===\n");
  int n = 9;
  double* M = test_matrix(n);
  double q[9] = {-1., 0.8, -0.3, -0.5, -0.6, 0.4, 0.7, 0.2, -0.1};
  double mu[3] = {0.4, 0.6, 0.3};
  BlockSpec one[1] = {{MIXED_MODEL_FRICTION_3D, 0, 3, 3, mu, NULL}};
  BlockSpec per_contact[3];
  for (int c = 0; c < 3; c++)
    per_contact[c] = (BlockSpec){MIXED_MODEL_FRICTION_3D, 3 * c, 3, 1, &mu[c], NULL};
  int failures = solve_and_compare("friction3d_contacts", n, M, q, 3, per_contact, NM_DENSE, 1,
                                   one, NM_DENSE);
  printf("=== Friction 3D, 3 contacts: dense vs sparse M ===\n");
  failures +=
      solve_and_compare("friction3d_sparse", n, M, q, 1, one, NM_DENSE, 1, one, NM_SPARSE);
  free(M);
  return failures;
}

/* Blocks of different dimensions: Friction 3D (1 contact), Friction 2D (2 contacts),
   Mohr-Coulomb (1 cone). The element rows are not multiples of the element dimension, so the
   generic NM_extract_diag_block is used. Dense vs sparse M. */
static int test_heterogeneous_blocks(void) {
  printf("=== Friction 3D + Friction 2D (2 contacts) + Mohr-Coulomb: dense vs sparse M ===\n");
  int n = 10;
  double* M = test_matrix(n);
  double q[10] = {-1., 0.6, -0.2, -0.9, 0.5, -1.1, -0.3, -1., -0.8, 0.3};
  double mu3[1] = {0.4}, mu2[2] = {0.3, 0.5}, eta[1] = {0.5}, theta[1] = {0.2};
  BlockSpec blocks[3] = {{MIXED_MODEL_FRICTION_3D, 0, 3, 1, mu3, NULL},
                         {MIXED_MODEL_FRICTION_2D, 3, 2, 2, mu2, NULL},
                         {MIXED_MODEL_DRUCKER_PRAGER, 7, 3, 1, eta, theta}};
  int failures = solve_and_compare("heterogeneous_blocks", n, M, q, 3, blocks, NM_DENSE, 3,
                                   blocks, NM_SPARSE);
  free(M);
  return failures;
}

/* Rolling friction 2D contact (dimension 3) coupled with a friction 2D contact. */
static int test_rolling_friction2d(void) {
  printf("=== Rolling friction 2D + Friction 2D ===\n");
  int n = 5;
  double* M = test_matrix(n);
  double q[5] = {-1., 0.5, 0.3, -0.7, 0.4};
  double mu[1] = {0.4}, mur[1] = {0.1}, mu2[1] = {0.3};
  BlockSpec blocks[2] = {{MIXED_MODEL_ROLLING_FRICTION_2D, 0, 3, 1, mu, mur},
                         {MIXED_MODEL_FRICTION_2D, 3, 2, 1, mu2, NULL}};
  double z[5] = {0.}, w[5] = {0.};
  int info = solve_mixed(n, M, q, 2, blocks, z, w, NM_DENSE);
  int failures =
      check_solution("rolling_friction2d", info, n, M, q, 2, blocks, z, w, NULL, 0.);
  free(M);
  return failures;
}

/* Rolling friction 3D contact (dimension 5) coupled with a friction 3D contact. */
static int test_rolling_friction3d(void) {
  printf("=== Rolling friction 3D + Friction 3D ===\n");
  int n = 8;
  double* M = test_matrix(n);
  double q[8] = {-1., 0.5, -0.2, 0.3, -0.1, -0.8, 0.4, 0.2};
  double mu[1] = {0.4}, mur[1] = {0.1}, mu3[1] = {0.3};
  BlockSpec blocks[2] = {{MIXED_MODEL_ROLLING_FRICTION_3D, 0, 5, 1, mu, mur},
                         {MIXED_MODEL_FRICTION_3D, 5, 3, 1, mu3, NULL}};
  double z[8] = {0.}, w[8] = {0.};
  int info = solve_mixed(n, M, q, 2, blocks, z, w, NM_DENSE);
  int failures =
      check_solution("rolling_friction3d", info, n, M, q, 2, blocks, z, w, NULL, 0.);
  free(M);
  return failures;
}

int main(void) {
  int failures = 0;
  failures += test_friction2d_contacts();
  failures += test_friction3d_contacts();
  failures += test_heterogeneous_blocks();
  failures += test_rolling_friction2d();
  failures += test_rolling_friction3d();
  printf("%d test(s) FAILED\n", failures);
  return failures ? 1 : 0;
}
