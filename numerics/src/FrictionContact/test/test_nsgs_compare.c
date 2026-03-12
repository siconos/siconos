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

/*!\file test_nsgs_compare.c
 * \brief Detailed comparison test between fc3d_nsgs and fc3d_nsgs_generic
 *
 * This test runs both solvers on a small problem and compares their
 * iteration-by-iteration behavior using stats callbacks.
 */

#define _POSIX_C_SOURCE 200809L  /* for clock_gettime and timespec */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "NumericsFwd.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#include "numerics_verbose.h"

/* Maximum number of iterations to track */
#define MAX_ITER_TRACK 1000

/* Structure to hold iteration data */
typedef struct {
  int iter;
  double error;
  double reaction_norm;
} IterationData;

/* Global arrays to store iteration history */
static IterationData original_iters[MAX_ITER_TRACK];
static IterationData generic_iters[MAX_ITER_TRACK];
static int original_count = 0;
static int generic_count = 0;

/* Stats callback for original fc3d_nsgs */
static void stats_callback_original(void* problem, SolverOptions* options, 
                                    double* reaction, double* velocity, double error) {
  (void)problem;
  (void)options;
  (void)velocity;
  
  if (original_count < MAX_ITER_TRACK) {
    original_iters[original_count].iter = options->iparam[SICONOS_IPARAM_ITER_DONE];
    original_iters[original_count].error = error;
    
    /* Compute reaction norm */
    int nc = ((FrictionContactProblem*)problem)->numberOfContacts;
    double rnorm = 0.0;
    for (int i = 0; i < nc * 3; i++) {
      rnorm += reaction[i] * reaction[i];
    }
    original_iters[original_count].reaction_norm = sqrt(rnorm);
    original_count++;
  }
}

/* Stats callback for generic fc3d_nsgs_generic */
static void stats_callback_generic(void* problem, SolverOptions* options,
                                   double* reaction, double* velocity, double error) {
  (void)problem;
  (void)options;
  (void)velocity;
  
  if (generic_count < MAX_ITER_TRACK) {
    generic_iters[generic_count].iter = options->iparam[SICONOS_IPARAM_ITER_DONE];
    generic_iters[generic_count].error = error;
    
    /* Compute reaction norm */
    int nc = ((FrictionContactProblem*)problem)->numberOfContacts;
    double rnorm = 0.0;
    for (int i = 0; i < nc * 3; i++) {
      rnorm += reaction[i] * reaction[i];
    }
    generic_iters[generic_count].reaction_norm = sqrt(rnorm);
    generic_count++;
  }
}

/* Create a simple test problem with 5 contacts */
static FrictionContactProblem* create_test_problem(void) {
  int nc = 5;
  int n = nc * 3;
  
  FrictionContactProblem* problem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
  problem->numberOfContacts = nc;
  problem->dimension = 3;
  problem->q = (double*)malloc(n * sizeof(double));
  problem->mu = (double*)malloc(nc * sizeof(double));
  
  /* Create a simple dense matrix with diagonal dominance */
  double* M_data = (double*)malloc(n * n * sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        M_data[i * n + j] = 2.0;  /* Diagonal */
      } else if (i/3 == j/3) {
        M_data[i * n + j] = 0.3;  /* Same contact coupling */
      } else {
        M_data[i * n + j] = 0.05;  /* Weak coupling between contacts */
      }
    }
  }
  
  problem->M = NM_create_from_data(NM_DENSE, n, n, M_data);
  
  /* Set q vector */
  for (int i = 0; i < n; i++) {
    problem->q[i] = -1.0 + 0.1 * (i % 3);
  }
  
  /* Set friction coefficients */
  for (int i = 0; i < nc; i++) {
    problem->mu[i] = 0.5;
  }
  
  return problem;
}

/* Print iteration comparison */
static void print_iteration_comparison(void) {
  printf("\n");
  printf("=================================================================\n");
  printf("Iteration-by-iteration Comparison\n");
  printf("=================================================================\n");
  printf("%-6s | %-20s | %-20s | %-12s\n", 
         "Iter", "Original Error", "Generic Error", "Diff");
  printf("-------+----------------------+----------------------+-------------\n");
  
  int max_iter = (original_count > generic_count) ? original_count : generic_count;
  
  for (int i = 0; i < max_iter && i < 50; i++) {  /* Limit to 50 iterations for display */
    double orig_err = (i < original_count) ? original_iters[i].error : 0.0;
    double gen_err = (i < generic_count) ? generic_iters[i].error : 0.0;
    double diff = fabs(orig_err - gen_err);
    
    char orig_str[32], gen_str[32];
    if (i < original_count) {
      snprintf(orig_str, 32, "%.6e", orig_err);
    } else {
      snprintf(orig_str, 32, "%s", "N/A");
    }
    
    if (i < generic_count) {
      snprintf(gen_str, 32, "%.6e", gen_err);
    } else {
      snprintf(gen_str, 32, "%s", "N/A");
    }
    
    printf("%-6d | %-20s | %-20s | %.6e\n", i+1, orig_str, gen_str, diff);
  }
  
  if (max_iter > 50) {
    printf("... (truncated after 50 iterations, total: %d)\n", max_iter);
  }
  
  printf("=================================================================\n");
  printf("Total iterations: original=%d, generic=%d\n", original_count, generic_count);
  printf("=================================================================\n");
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  
  printf("============================================================\n");
  printf("NSGS Detailed Comparison Test\n");
  printf("============================================================\n");
  
  /* Create test problem */
  FrictionContactProblem* problem = create_test_problem();
  int nc = problem->numberOfContacts;
  
  printf("Problem size: %d contacts (%d variables)\n", nc, nc * 3);
  printf("\n");
  
  /* Setup solver options */
  SolverOptions* options_orig = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options_orig->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options_orig->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  
  /* Setup internal solver (projection on cone) */
  solver_options_update_internal(options_orig, 0, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  options_orig->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options_orig->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  
  /* Copy options for generic solver */
  SolverOptions* options_gen = solver_options_copy(options_orig);
  
  /* Allocate solution arrays */
  double* r_orig = (double*)calloc(nc * 3, sizeof(double));
  double* v_orig = (double*)calloc(nc * 3, sizeof(double));
  double* r_gen = (double*)calloc(nc * 3, sizeof(double));
  double* v_gen = (double*)calloc(nc * 3, sizeof(double));
  
  /* Run original fc3d_nsgs */
  printf("Running original fc3d_nsgs...\n");
  original_count = 0;
  int info_orig;
  fc3d_nsgs(problem, r_orig, v_orig, &info_orig, options_orig);
  
  /* Run generic fc3d_nsgs_generic */
  printf("Running generic fc3d_nsgs_generic...\n");
  generic_count = 0;
  int info_gen;
  fc3d_nsgs_generic(problem, r_gen, v_gen, &info_gen, options_gen);
  
  /* Print iteration comparison */
  print_iteration_comparison();
  
  /* Compare final solutions */
  printf("\n");
  printf("============================================================\n");
  printf("Final Solution Comparison\n");
  printf("============================================================\n");
  
  double diff_r = 0.0, diff_v = 0.0;
  double norm_r_orig = 0.0, norm_r_gen = 0.0;
  
  for (int i = 0; i < nc * 3; i++) {
    diff_r += (r_orig[i] - r_gen[i]) * (r_orig[i] - r_gen[i]);
    diff_v += (v_orig[i] - v_gen[i]) * (v_orig[i] - v_gen[i]);
    norm_r_orig += r_orig[i] * r_orig[i];
    norm_r_gen += r_gen[i] * r_gen[i];
  }
  diff_r = sqrt(diff_r);
  diff_v = sqrt(diff_v);
  norm_r_orig = sqrt(norm_r_orig);
  norm_r_gen = sqrt(norm_r_gen);
  
  printf("Final reaction norm (original): %.6e\n", norm_r_orig);
  printf("Final reaction norm (generic):  %.6e\n", norm_r_gen);
  printf("Difference in reaction:         %.6e\n", diff_r);
  printf("Difference in velocity:         %.6e\n", diff_v);
  printf("Convergence: original=%d, generic=%d\n", info_orig, info_gen);
  printf("Iterations: original=%d, generic=%d\n", 
         options_orig->iparam[SICONOS_IPARAM_ITER_DONE],
         options_gen->iparam[SICONOS_IPARAM_ITER_DONE]);
  
  /* Print first few reaction values */
  printf("\nFirst 6 reaction values:\n");
  printf("%-8s | %-15s | %-15s | %-12s\n", "Index", "Original", "Generic", "Diff");
  printf("---------+-----------------+-----------------+-------------\n");
  for (int i = 0; i < 6 && i < nc * 3; i++) {
    printf("%-8d | %15.6e | %15.6e | %12.6e\n", 
           i, r_orig[i], r_gen[i], fabs(r_orig[i] - r_gen[i]));
  }
  
  /* Performance comparison */
  printf("\n");
  printf("============================================================\n");
  printf("Performance Comparison (100 runs)\n");
  printf("============================================================\n");
  
  /* Time original */
  struct timespec start, end;
  double time_orig = 0.0;
  for (int run = 0; run < 100; run++) {
    double* r_tmp = (double*)calloc(nc * 3, sizeof(double));
    double* v_tmp = (double*)calloc(nc * 3, sizeof(double));
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    fc3d_nsgs(problem, r_tmp, v_tmp, &info_orig, options_orig);
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    time_orig += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    
    free(r_tmp);
    free(v_tmp);
  }
  time_orig /= 100.0;
  
  /* Time generic */
  double time_generic = 0.0;
  for (int run = 0; run < 100; run++) {
    double* r_tmp = (double*)calloc(nc * 3, sizeof(double));
    double* v_tmp = (double*)calloc(nc * 3, sizeof(double));
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    fc3d_nsgs_generic(problem, r_tmp, v_tmp, &info_gen, options_gen);
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    time_generic += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    
    free(r_tmp);
    free(v_tmp);
  }
  time_generic /= 100.0;
  
  printf("Average time (original): %.6f ms\n", time_orig * 1000);
  printf("Average time (generic):  %.6f ms\n", time_generic * 1000);
  printf("Overhead: %.1f%%\n", ((time_generic - time_orig) / time_orig) * 100.0);
  
  printf("============================================================\n");
  
  /* Cleanup */
  free(r_orig);
  free(v_orig);
  free(r_gen);
  free(v_gen);
  solver_options_delete(options_orig);
  solver_options_delete(options_gen);
  frictionContactProblem_free(problem);
  
  return (diff_r < 1e-6 && diff_v < 1e-6) ? 0 : 1;
}
