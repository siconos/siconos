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

/*!
 * \file test_cohesive_friction_3d_simple.c
 * \brief Simple test for cohesive friction 3D solver
 *
 * This test creates a simple problem with:
 * - 1 contact point
 * - 1 cohesive point
 * - All matrices (W, V, U, X) equal to identity
 * - Simple q_v and q_u vectors
 *
 * The expected behavior is that the solver should converge to a solution
 * where the reaction satisfies the friction cone condition.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CohesiveFrictionContactProblem.h"
#include "CohesiveFrictionContact_options.h"
#include "NonSmoothGaussSeidel_options.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SolverOptions.h"
#include "cohesive_friction_3d_driver.h"
#include "numerics_verbose.h"

/* Helper function to create an identity matrix */
static NumericsMatrix* create_identity_matrix(int size) {
  NumericsMatrix* M = NM_create(NM_SPARSE, size, size);
  NM_triplet_alloc(M, 0);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      NM_zentry(M, i, j, (i == j) ? 1.0 : 0.0, 1e-12);
    }
  }
  return M;
}

/* Helper function to print a vector */
static void print_vector(const char* name, double* vec, int size) {
  printf("%s = [", name);
  for (int i = 0; i < size; i++) {
    printf(" %10.6e", vec[i]);
    if (i < size - 1) printf(",");
  }
  printf(" ]\n");
}
static CohesiveFrictionContactProblem* build_problem1x1() {
  /* Test parameters */
  int dim = 3;          // 3D problem
  int nc = 1;           // 1 contact point
  int n_coh = 1;        // 1 cohesive point
  int m = dim * nc;     // size of velocity block = 3
  int n = dim * n_coh;  // size of cohesive block = 3
  int M_size = m + n;   // total size = 6

  printf("Problem setup:\n");
  printf("  Dimension: %d (3D)\n", dim);
  printf("  Number of contacts: %d\n", nc);
  printf("  Number of cohesive points: %d\n", n_coh);
  printf("  Total problem size: %d\n\n", M_size);

  /* Create problem structure */
  CohesiveFrictionContactProblem* problem = cohesiveFrictionContactProblem_new();
  problem->dimension = dim;
  problem->numberOfContacts = nc;
  problem->numberOfCohesivePoints = n_coh;

  /* Create block matrices (all identity) */
  printf("Creating block matrices (all identity)...\n");

  // W: contact-contact block (3x3)
  problem->W = create_identity_matrix(m);
  printf("  W: %dx%d identity matrix\n", m, m);

  // V: contact-cohesive coupling block (3x3)
  problem->V = create_identity_matrix(m);
  printf("  V: %dx%d identity matrix\n", m, n);

  // U: cohesive-contact coupling block (3x3)
  problem->U = create_identity_matrix(n);
  printf("  U: %dx%d identity matrix\n", n, m);

  // X: cohesive-cohesive block (3x3)
  problem->X = create_identity_matrix(n);
  printf("  X: %dx%d identity matrix\n", n, n);

  /* Create q vectors */
  printf("\nCreating q vectors...\n");

  // q_v: contact part (initial gap/velocity)
  problem->q_v = (double*)malloc(m * sizeof(double));
  problem->q_v[0] = 0.1;  // normal gap (negative = penetration)
  problem->q_v[1] = 1.0;  // tangent 1
  problem->q_v[2] = 1.0;  // tangent 2
  print_vector("  q_v", problem->q_v, m);

  // q_u: cohesive part (initial displacement)
  problem->q_u = (double*)malloc(n * sizeof(double));
  double h = 1;
  double u_k = 0.5;
  problem->q_u[0] = u_k + h * problem->q_v[0];  // normal displacement
  problem->q_u[1] = 0.0;                        // tangent 1
  problem->q_u[2] = 0.0;                        // tangent 2
  print_vector("  q_u", problem->q_u, n);

  /* Friction coefficient */
  problem->mu = (double*)malloc(nc * sizeof(double));
  problem->mu[0] = 0.5;  // friction coefficient
  printf("\nFriction coefficient: mu = %.2f\n", problem->mu[0]);

  /* Cohesion intensities */
  problem->c_n = (double*)malloc(n_coh * sizeof(double));
  problem->c_t = (double*)malloc(n_coh * sizeof(double));
  problem->c_n[0] = -1.0;  // normal cohesion
  problem->c_t[0] = 0.0;   // tangent cohesion
  printf("Normal cohesion: c_n = %.2f\n", problem->c_n[0]);
  printf("Tangent cohesion: c_t = %.2f\n", problem->c_t[0]);

  /* Build global M and q from blocks */
  printf("\nBuilding global M and q from blocks...\n");
  int build_status = cohesiveFrictionContactProblem_build_M_q_from_blocks(problem);
  if (build_status != 0) {
    fprintf(stderr, "ERROR: Failed to build M and q from blocks\n");
    return NULL;
  }
  printf("  M: %dx%d block matrix [[W,V],[U,X]]\n", M_size, M_size);
  printf("  q: %d vector [q_v, q_u]\n", M_size);

  return problem;
}

static CohesiveFrictionContactProblem* build_problem0x1() {
  /* Test parameters */
  int dim = 3;          // 3D problem
  int nc = 0;           // 1 contact point
  int n_coh = 1;        // 1 cohesive point
  int m = dim * nc;     // size of velocity block = 3
  int n = dim * n_coh;  // size of cohesive block = 3
  int M_size = m + n;   // total size = 6

  printf("Problem setup:\n");
  printf("  Dimension: %d (3D)\n", dim);
  printf("  Number of contacts: %d\n", nc);
  printf("  Number of cohesive points: %d\n", n_coh);
  printf("  Total problem size: %d\n\n", M_size);

  /* Create problem structure */
  CohesiveFrictionContactProblem* problem = cohesiveFrictionContactProblem_new();
  problem->dimension = dim;
  problem->numberOfContacts = nc;
  problem->numberOfCohesivePoints = n_coh;

  /* Create block matrices (all identity) */
  printf("Creating block matrices (all identity)...\n");

  /* // W: contact-contact block (3x3) */
  /* problem->W = create_identity_matrix(m); */
  /* printf("  W: %dx%d identity matrix\n", m, m); */

  /* // V: contact-cohesive coupling block (3x3) */
  /* problem->V = create_identity_matrix(m); */
  /* printf("  V: %dx%d identity matrix\n", m, n); */

  /* // U: cohesive-contact coupling block (3x3) */
  /* problem->U = create_identity_matrix(n); */
  /* printf("  U: %dx%d identity matrix\n", n, m); */

  // X: cohesive-cohesive block (3x3)
  problem->X = create_identity_matrix(n);
  printf("  X: %dx%d identity matrix\n", n, n);

  /* Create q vectors */
  printf("\nCreating q vectors...\n");

  /* // q_v: contact part (initial gap/velocity) */
  /* problem->q_v = (double*)malloc(m * sizeof(double)); */
  /* problem->q_v[0] = 0.1;  // normal gap (negative = penetration) */
  /* problem->q_v[1] = 1.0;   // tangent 1 */
  /* problem->q_v[2] = 1.0;   // tangent 2 */
  /* print_vector("  q_v", problem->q_v, m); */

  // q_u: cohesive part (initial displacement)
  problem->q_u = (double*)malloc(n * sizeof(double));
  double h = 1;
  double u_k = 0.5;
  problem->q_u[0] = u_k;  // normal displacement
  problem->q_u[1] = 0.0;  // tangent 1
  problem->q_u[2] = 0.0;  // tangent 2
  print_vector("  q_u", problem->q_u, n);

  /* Friction coefficient */
  /* problem->mu = (double*)malloc(nc * sizeof(double)); */
  /* problem->mu[0] = 0.5;    // friction coefficient */
  /* printf("\nFriction coefficient: mu = %.2f\n", problem->mu[0]); */

  /* Cohesion intensities */
  problem->c_n = (double*)malloc(n_coh * sizeof(double));
  problem->c_t = (double*)malloc(n_coh * sizeof(double));
  problem->c_n[0] = -1.0;  // normal cohesion
  problem->c_t[0] = 0.0;   // tangent cohesion
  printf("Normal cohesion: c_n = %.2f\n", problem->c_n[0]);
  printf("Tangent cohesion: c_t = %.2f\n", problem->c_t[0]);

  /* Build global M and q from blocks */
  printf("\nBuilding global M and q from blocks...\n");
  int build_status = cohesiveFrictionContactProblem_build_M_q_from_blocks(problem);
  if (build_status != 0) {
    fprintf(stderr, "ERROR: Failed to build M and q from blocks\n");
    return NULL;
  }
  printf("  M: %dx%d block matrix [[W,V],[U,X]]\n", M_size, M_size);
  printf("  q: %d vector [q_v, q_u]\n", M_size);

  return problem;
}


static CohesiveFrictionContactProblem* build_problem1x0() {
  /* Test parameters */
  int dim = 3;          // 3D problem
  int nc = 1;           // 1 contact point
  int n_coh = 0;        // 1 cohesive point
  int m = dim * nc;     // size of velocity block = 3
  int n = dim * n_coh;  // size of cohesive block = 3
  int M_size = m + n;   // total size = 6

  printf("Problem setup:\n");
  printf("  Dimension: %d (3D)\n", dim);
  printf("  Number of contacts: %d\n", nc);
  printf("  Number of cohesive points: %d\n", n_coh);
  printf("  Total problem size: %d\n\n", M_size);

  /* Create problem structure */
  CohesiveFrictionContactProblem* problem = cohesiveFrictionContactProblem_new();
  problem->dimension = dim;
  problem->numberOfContacts = nc;
  problem->numberOfCohesivePoints = n_coh;

  /* Create block matrices (all identity) */
  printf("Creating block matrices (all identity)...\n");

  // W: contact-contact block (3x3)
  problem->W = create_identity_matrix(m);
  printf("  W: %dx%d identity matrix\n", m, m);

  /* // V: contact-cohesive coupling block (3x3) */
  /* problem->V = create_identity_matrix(m); */
  /* printf("  V: %dx%d identity matrix\n", m, n); */

  /* // U: cohesive-contact coupling block (3x3) */
  /* problem->U = create_identity_matrix(n); */
  /* printf("  U: %dx%d identity matrix\n", n, m); */

  /* // X: cohesive-cohesive block (3x3) */
  /* problem->X = create_identity_matrix(n); */
  /* printf("  X: %dx%d identity matrix\n", n, n); */

  /* Create q vectors */
  printf("\nCreating q vectors...\n");

  // q_v: contact part (initial gap/velocity)
  problem->q_v = (double*)malloc(m * sizeof(double));
  problem->q_v[0] = -0.1;  // normal gap (negative = penetration)
  problem->q_v[1] = 1.0;  // tangent 1
  problem->q_v[2] = 1.0;  // tangent 2
  print_vector("  q_v", problem->q_v, m);

  /* // q_u: cohesive part (initial displacement) */
  /* problem->q_u = (double*)malloc(n * sizeof(double)); */
  /* double h = 1; */
  /* double u_k = 0.5; */
  /* problem->q_u[0] = u_k + h * problem->q_v[0];  // normal displacement */
  /* problem->q_u[1] = 0.0;                        // tangent 1 */
  /* problem->q_u[2] = 0.0;                        // tangent 2 */
  /* print_vector("  q_u", problem->q_u, n); */

  /* Friction coefficient */
  problem->mu = (double*)malloc(nc * sizeof(double));
  problem->mu[0] = 0.5;  // friction coefficient
  printf("\nFriction coefficient: mu = %.2f\n", problem->mu[0]);

  /* /\* Cohesion intensities *\/ */
  /* problem->c_n = (double*)malloc(n_coh * sizeof(double)); */
  /* problem->c_t = (double*)malloc(n_coh * sizeof(double)); */
  /* problem->c_n[0] = -1.0;  // normal cohesion */
  /* problem->c_t[0] = 0.0;   // tangent cohesion */
  /* printf("Normal cohesion: c_n = %.2f\n", problem->c_n[0]); */
  /* printf("Tangent cohesion: c_t = %.2f\n", problem->c_t[0]); */

  /* Build global M and q from blocks */
  printf("\nBuilding global M and q from blocks...\n");
  int build_status = cohesiveFrictionContactProblem_build_M_q_from_blocks(problem);
  if (build_status != 0) {
    fprintf(stderr, "ERROR: Failed to build M and q from blocks\n");
    return NULL;
  }
  printf("  M: %dx%d block matrix [[W,V],[U,X]]\n", M_size, M_size);
  printf("  q: %d vector [q_v, q_u]\n", M_size);

  return problem;
}

int test_problem(CohesiveFrictionContactProblem* problem) {
  int M_size =
      (problem->numberOfContacts + problem->numberOfCohesivePoints) * problem->dimension;
  int nc = problem->numberOfContacts;
  int ncoh = problem->numberOfCohesivePoints;

  /* Display problem */
  printf("\n");
  cohesiveFrictionContact_display(problem);

  /* Allocate solution vectors */
  double* reaction = (double*)calloc(M_size, sizeof(double));
  double* velocity = (double*)calloc(M_size, sizeof(double));

  /* Create solver options */
  printf("\nSetting up solver...\n");
  SolverOptions* options = solver_options_create(SICONOS_COHESIVE_FRICTION_3D_NSGS);
  if (!options) {
    fprintf(stderr, "ERROR: Failed to create solver options\n");
    return 1;
  }

  /* Set solver parameters */
  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->iparam[SICONOS_NSGS_FREEZING_CONTACT] = 0;

  printf("  Solver: NSGS (Non-Smooth Gauss-Seidel)\n");
  printf("  Tolerance: %.2e\n", options->dparam[SICONOS_DPARAM_TOL]);
  printf("  Max iterations: %d\n", options->iparam[SICONOS_IPARAM_MAX_ITER]);

  /* Solve the problem */
  printf("\nSolving...\n");
  int info = cohesive_friction_3d_driver(problem, reaction, velocity, options);

  /* Print results */
  printf("\n=================================================================\n");
  printf("=== Results =====================================================\n");
  printf("=================================================================\n");

  if (info == 0) {
    printf("Solver converged successfully!\n\n");
    for (int i = 0; i < nc; i++) {
      printf("Contact reaction (p):\n");
      printf("  r_n  = %12.6e (normal)\n", reaction[0 + 3 * i]);
      printf("  r_t1 = %12.6e (tangent 1)\n", reaction[1 + 3 * i]);
      printf("  r_t2 = %12.6e (tangent 2)\n", reaction[2 + 3 * i]);
      printf("\nContact velocity (v):\n");
      printf("  v_n  = %12.6e (normal)\n", velocity[0 + 3 * i]);
      printf("  v_t1 = %12.6e (tangent 1)\n", velocity[1 + 3 * i]);
      printf("  v_t2 = %12.6e (tangent 2)\n", velocity[2 + 3 * i]);
      /* Verify friction cone condition */
      double r_n = reaction[0];
      double r_t_norm = sqrt(reaction[1] * reaction[1] + reaction[2] * reaction[2]);
      double mu_r_n = problem->mu[0] * r_n;

      printf("\nVerification:\n");
      printf("  Normal reaction r_n = %.6e\n", r_n);
      printf("  Tangential reaction ||r_t|| = %.6e\n", r_t_norm);
      printf("  Friction bound mu*r_n = %.6e\n", mu_r_n);
      printf("  Friction cone condition: ||r_t|| <= mu*r_n : %s\n",
             (r_t_norm <= mu_r_n + 1e-10) ? "SATISFIED" : "VIOLATED");

      /* Check complementarity */
      double u_n = velocity[0];
      double comp = r_n * u_n;
      printf("  Complementarity r_n * u_n = %.6e (should be ~0)\n", comp);
    }
    for (int k = nc; k < nc + ncoh; k++) {
      int i = k - nc;
      printf("\nCohesive reaction (r_coh):\n");
      printf("  r_n  = %12.6e (normal)\n", reaction[0 + 3 * i]);
      printf("  r_t1 = %12.6e (tangent 1)\n", reaction[1 + 3 * i]);
      printf("  r_t2 = %12.6e (tangent 2)\n", reaction[2 + 3 * i]);

      printf("\nCohesive displacement rate (u):\n");
      printf("  u_n  = %12.6e (normal)\n", velocity[0 + 3 * i]);
      printf("  u_t1 = %12.6e (tangent 1)\n", velocity[1 + 3 * i]);
      printf("  u_t2 = %12.6e (tangent 2)\n", velocity[2 + 3 * i]);
    }

 

  } else {
    printf("Solver failed with error code: %d\n", info);
  }

  printf("\nIterations: %d\n", options->iparam[SICONOS_IPARAM_ITER_DONE]);
  printf("Final residual: %.6e\n", options->dparam[SICONOS_DPARAM_RESIDU]);

  /* Cleanup */
  solver_options_delete(options);

  free(reaction);
  free(velocity);
  return info;
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  printf("=================================================================\n");
  printf("=== Cohesive Friction 3D Simple Test ============================\n");
  printf("=================================================================\n\n");
  CohesiveFrictionContactProblem* problem = build_problem1x1();
  int info = test_problem(problem);
  cohesiveFrictionContactProblem_free(problem);

  CohesiveFrictionContactProblem* problem_0x1 = build_problem0x1();
  info += test_problem(problem_0x1);
  cohesiveFrictionContactProblem_free(problem_0x1);
  
  CohesiveFrictionContactProblem* problem_1x0 = build_problem1x0();
  info += test_problem(problem_1x0);
  cohesiveFrictionContactProblem_free(problem_1x0);

  printf("\n=================================================================\n");
  printf("=== Test %s ======================================\n",
         (info == 0) ? "PASSED" : "FAILED");
  printf("=================================================================\n");

  return info;
}
