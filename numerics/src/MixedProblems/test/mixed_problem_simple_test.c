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
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MixedProblem.h"
#include "MixedProblem_NSGS.h"
#include "MixedProblems_options.h"
#include "NumericsMatrix.h"
#include "NumericsVerbose.h"
#include "SolverOptions.h"
/*
 * Simple test for MixedProblem with:
 * - M = Identity matrix
 * - q = vector of -1
 * - One block of each supported model type
 */

static int test_friction2d_only(void) {
  printf("\n=== Test: Friction 2D only ===\n");

  int dim = 2;
  int ne = 1;
  int total_dim = dim * ne;

  /* Create M = Identity */
  NumericsMatrix* M = NM_create(NM_DENSE, total_dim, total_dim);
  for (int i = 0; i < total_dim; i++) {
    for (int j = 0; j < total_dim; j++) {
      M->dense[i * total_dim + j] = (i == j) ? 1.0 : 0.0;
    }
  }

  /* Create q = -1 */
  double* q = (double*)calloc(total_dim, sizeof(double));
  for (int i = 0; i < total_dim; i++) {
    q[i] = -1.0;
  }

  /* Create MixedProblem with 1 block */
  MixedProblem* problem = mixedProblem_new_with_data(1, M, q);

  /* Set block: Friction 2D */
  mixedProblem_set_block(problem, 0, 0, 2, 1, MIXED_MODEL_FRICTION_2D);
  double mu = 0.3;
  mixedProblem_set_friction2d_params(problem, 0, &mu);

  /* Display problem */
  mixedProblem_display(problem);

  /* Solve with NSGS */
  SolverOptions options;
  mixedProblem_NSGS_set_default_solver_options(problem, &options);
  options.iparam[SICONOS_MIXED_IPARAM_MAX_ITER] = 100;
  options.dparam[SICONOS_MIXED_DPARAM_TOL] = 1e-6;

  double* z = (double*)calloc(total_dim, sizeof(double));
  double* w = (double*)calloc(total_dim, sizeof(double));
  int info;

  mixedProblem_driver(problem, z, w, &options);
  info = options.dparam[SICONOS_MIXED_DPARAM_RESIDU] < options.dparam[SICONOS_MIXED_DPARAM_TOL]
             ? 0
             : 1;

  printf("\nSolution:\n");
  printf("  z = [%12.8e, %12.8e]\n", z[0], z[1]);
  printf("  w = [%12.8e, %12.8e]\n", w[0], w[1]);
  printf("  iterations = %d\n", options.iparam[SICONOS_MIXED_IPARAM_ITER_DONE]);
  printf("  residual = %12.8e\n", options.dparam[SICONOS_MIXED_DPARAM_RESIDU]);

  /* Cleanup */
  free(z);
  free(w);
  NM_free(M);
  free(q);
  mixedProblem_free(problem);
  solver_options_delete(&options);

  printf("Test %s\n", info == 0 ? "PASSED" : "FAILED");
  return info;
}

static int test_vonmises_only(void) {
  printf("\n=== Test: Von Mises only ===\n");

  int dim = 3;
  int ne = 1;
  int total_dim = dim * ne;

  /* Create M = Identity */
  NumericsMatrix* M = NM_create(NM_DENSE, total_dim, total_dim);
  for (int i = 0; i < total_dim; i++) {
    for (int j = 0; j < total_dim; j++) {
      M->dense[i * total_dim + j] = (i == j) ? 1.0 : 0.0;
    }
  }

  /* Create q = -1 */
  double* q = (double*)calloc(total_dim, sizeof(double));
  for (int i = 0; i < total_dim; i++) {
    q[i] = -1.0;
  }

  /* Create MixedProblem with 1 block */
  MixedProblem* problem = mixedProblem_new_with_data(1, M, q);

  /* Set block: Von Mises */
  mixedProblem_set_block(problem, 0, 0, 3, 1, MIXED_MODEL_VON_MISES);
  double sigma_y = 100.0;
  mixedProblem_set_vonmises_params(problem, 0, &sigma_y);

  /* Display problem */
  mixedProblem_display(problem);

  /* Solve with NSGS */
  SolverOptions options;
  mixedProblem_NSGS_set_default_solver_options(problem, &options);
  options.iparam[SICONOS_MIXED_IPARAM_MAX_ITER] = 100;
  options.dparam[SICONOS_MIXED_DPARAM_TOL] = 1e-6;

  double* z = (double*)calloc(total_dim, sizeof(double));
  double* w = (double*)calloc(total_dim, sizeof(double));
  int info;

  mixedProblem_driver(problem, z, w, &options);
  info = options.dparam[SICONOS_MIXED_DPARAM_RESIDU] < options.dparam[SICONOS_MIXED_DPARAM_TOL]
             ? 0
             : 1;

  printf("\nSolution:\n");
  printf("  z = [%12.8e, %12.8e, %12.8e]\n", z[0], z[1], z[2]);
  printf("  w = [%12.8e, %12.8e, %12.8e]\n", w[0], w[1], w[2]);
  printf("  iterations = %d\n", options.iparam[SICONOS_MIXED_IPARAM_ITER_DONE]);
  printf("  residual = %12.8e\n", options.dparam[SICONOS_MIXED_DPARAM_RESIDU]);

  /* Cleanup */
  free(z);
  free(w);
  NM_free(M);
  free(q);
  mixedProblem_free(problem);
  solver_options_delete(&options);

  printf("Test %s\n", info == 0 ? "PASSED" : "FAILED");
  return info;
}

static int test_friction2d_and_vonmises(void) {
  printf("\n=== Test: Friction 2D + Von Mises ===\n");

  /* Block 0: Friction 2D (dim=2)
   * Block 1: Von Mises (dim=3)
   * Total dimension = 5
   */
  int total_dim = 5;

  /* Create M = Identity */
  NumericsMatrix* M = NM_create(NM_DENSE, total_dim, total_dim);
  for (int i = 0; i < total_dim; i++) {
    for (int j = 0; j < total_dim; j++) {
      M->dense[i * total_dim + j] = (i == j) ? 1.0 : 0.0;
    }
  }

  /* Create q = -1 */
  double* q = (double*)calloc(total_dim, sizeof(double));
  for (int i = 0; i < total_dim; i++) {
    q[i] = -1.0;
  }

  /* Create MixedProblem with 2 blocks */
  MixedProblem* problem = mixedProblem_new_with_data(2, M, q);

  /* Set block 0: Friction 2D */
  mixedProblem_set_block(problem, 0, 0, 2, 1, MIXED_MODEL_FRICTION_2D);
  double mu = 0.3;
  mixedProblem_set_friction2d_params(problem, 0, &mu);

  /* Set block 1: Von Mises */
  mixedProblem_set_block(problem, 1, 2, 3, 1, MIXED_MODEL_VON_MISES);
  double sigma_y = 100.0;
  mixedProblem_set_vonmises_params(problem, 1, &sigma_y);

  /* Display problem */
  mixedProblem_display(problem);

  /* Solve with NSGS */
  SolverOptions options;
  mixedProblem_NSGS_set_default_solver_options(problem, &options);
  options.iparam[SICONOS_MIXED_IPARAM_MAX_ITER] = 100;
  options.dparam[SICONOS_MIXED_DPARAM_TOL] = 1e-6;

  double* z = (double*)calloc(total_dim, sizeof(double));
  double* w = (double*)calloc(total_dim, sizeof(double));
  int info;

  mixedProblem_driver(problem, z, w, &options);
  info = options.dparam[SICONOS_MIXED_DPARAM_RESIDU] < options.dparam[SICONOS_MIXED_DPARAM_TOL]
             ? 0
             : 1;

  printf("\nSolution:\n");
  printf("  Block 0 (Friction 2D): z = [%12.8e, %12.8e]\n", z[0], z[1]);
  printf("  Block 1 (Von Mises):   z = [%12.8e, %12.8e, %12.8e]\n", z[2], z[3], z[4]);
  printf("  iterations = %d\n", options.iparam[SICONOS_MIXED_IPARAM_ITER_DONE]);
  printf("  residual = %12.8e\n", options.dparam[SICONOS_MIXED_DPARAM_RESIDU]);

  /* Cleanup */
  free(z);
  free(w);
  NM_free(M);
  free(q);
  mixedProblem_free(problem);
  solver_options_delete(&options);

  printf("Test %s\n", info == 0 ? "PASSED" : "FAILED");
  return info;
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  numerics_set_verbose(1);

  printf("==================================================\n");
  printf("MixedProblems - Simple Test Suite\n");
  printf("M = Identity, q = -1\n");
  printf("==================================================\n");

  int total_info = 0;

  /* Test 1: Friction 2D only */
  total_info += test_friction2d_only();

  /* Test 2: Von Mises only */
  total_info += test_vonmises_only();

  /* Test 3: Friction 2D + Von Mises */
  total_info += test_friction2d_and_vonmises();

  printf("\n==================================================\n");
  printf("All tests %s\n", total_info == 0 ? "PASSED" : "FAILED");
  printf("==================================================\n");

  return total_info;
}
