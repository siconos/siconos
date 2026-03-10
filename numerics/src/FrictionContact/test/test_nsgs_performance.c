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

/*!\file test_nsgs_performance.c
 * \brief Performance comparison test between fc3d_nsgs and fc3d_nsgs_generic
 *
 * This test compares CPU time, iteration counts, and convergence behavior
 * between the original fc3d_nsgs solver and the new generic NSGS implementation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "NumericsFwd.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#include "numerics_verbose.h"

/* Get CPU time in seconds (user + system) */
static double get_cpu_time(void) {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec * 1e-6 +
         (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec * 1e-6;
}

/* Get wall clock time in seconds */
static double get_wall_time(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Run solver once and collect detailed statistics */
typedef struct {
  double cpu_time;
  double wall_time;
  int iterations;
  int info;
  double final_error;
  double reaction_norm;
} SolverStats;

static SolverStats run_solver_once(FrictionContactProblem* problem, 
                                   double* reaction, double* velocity,
                                   SolverOptions* options, int use_original) {
  SolverStats stats = {0};
  int nc = problem->numberOfContacts;
  
  /* Allocate working arrays */
  double* r_tmp = (double*)malloc(nc * 3 * sizeof(double));
  double* v_tmp = (double*)malloc(nc * 3 * sizeof(double));
  memcpy(r_tmp, reaction, nc * 3 * sizeof(double));
  memcpy(v_tmp, velocity, nc * 3 * sizeof(double));
  
  /* Run solver and measure time */
  double cpu_start = get_cpu_time();
  double wall_start = get_wall_time();
  
  if (use_original) {
    fc3d_nsgs(problem, r_tmp, v_tmp, &stats.info, options);
  } else {
    fc3d_nsgs_generic(problem, r_tmp, v_tmp, &stats.info, options);
  }
  
  stats.cpu_time = get_cpu_time() - cpu_start;
  stats.wall_time = get_wall_time() - wall_start;
  stats.iterations = options->iparam[SICONOS_IPARAM_ITER_DONE];
  stats.final_error = options->dparam[SICONOS_DPARAM_RESIDU];
  
  /* Compute reaction norm */
  double rnorm = 0.0;
  for (int i = 0; i < nc * 3; i++) {
    rnorm += r_tmp[i] * r_tmp[i];
  }
  stats.reaction_norm = sqrt(rnorm);
  
  free(r_tmp);
  free(v_tmp);
  
  return stats;
}

/* Run multiple trials and compute statistics */
typedef struct {
  double cpu_time_avg;
  double cpu_time_std;
  double wall_time_avg;
  double wall_time_std;
  int iter_avg;
  int iter_min;
  int iter_max;
  int converged_count;
} SolverBenchmark;

static SolverBenchmark benchmark_solver(FrictionContactProblem* problem, 
                                        double* reaction, double* velocity,
                                        SolverOptions* options, 
                                        int n_trials, int use_original,
                                        const char* solver_name) {
  SolverBenchmark bench = {0};
  SolverStats* stats = (SolverStats*)malloc(n_trials * sizeof(SolverStats));
  
  printf("  Benchmarking %s (%d trials)...\n", solver_name, n_trials);
  
  for (int i = 0; i < n_trials; i++) {
    stats[i] = run_solver_once(problem, reaction, velocity, options, use_original);
    if (stats[i].info == 0) bench.converged_count++;
  }
  
  /* Compute averages */
  double cpu_sum = 0.0, cpu_sq_sum = 0.0;
  double wall_sum = 0.0, wall_sq_sum = 0.0;
  int iter_sum = 0;
  bench.iter_min = stats[0].iterations;
  bench.iter_max = stats[0].iterations;
  
  for (int i = 0; i < n_trials; i++) {
    cpu_sum += stats[i].cpu_time;
    cpu_sq_sum += stats[i].cpu_time * stats[i].cpu_time;
    wall_sum += stats[i].wall_time;
    wall_sq_sum += stats[i].wall_time * stats[i].wall_time;
    iter_sum += stats[i].iterations;
    if (stats[i].iterations < bench.iter_min) bench.iter_min = stats[i].iterations;
    if (stats[i].iterations > bench.iter_max) bench.iter_max = stats[i].iterations;
  }
  
  bench.cpu_time_avg = cpu_sum / n_trials;
  bench.wall_time_avg = wall_sum / n_trials;
  bench.iter_avg = iter_sum / n_trials;
  
  /* Compute standard deviations */
  if (n_trials > 1) {
    bench.cpu_time_std = sqrt((cpu_sq_sum - n_trials * bench.cpu_time_avg * bench.cpu_time_avg) / (n_trials - 1));
    bench.wall_time_std = sqrt((wall_sq_sum - n_trials * bench.wall_time_avg * bench.wall_time_avg) / (n_trials - 1));
  }
  
  free(stats);
  return bench;
}

/* Print benchmark results */
static void print_benchmark_results(const char* name, SolverBenchmark* bench, int n_trials) {
  printf("\n%s:\n", name);
  printf("  CPU time:   %.6f +/- %.6f ms\n", bench->cpu_time_avg * 1000, bench->cpu_time_std * 1000);
  printf("  Wall time:  %.6f +/- %.6f ms\n", bench->wall_time_avg * 1000, bench->wall_time_std * 1000);
  printf("  Iterations: %d (min: %d, max: %d)\n", bench->iter_avg, bench->iter_min, bench->iter_max);
  printf("  Converged:  %d/%d\n", bench->converged_count, n_trials);
}

int main(int argc, char** argv) {
  /* Default test parameters */
  int n_trials = 10000;
  const char* data_file = "./data/KaplasTower-i1061-4.hdf5.dat";
  
  /* Parse command line arguments */
  if (argc > 1) {
    data_file = argv[1];
  }
  if (argc > 2) {
    n_trials = atoi(argv[2]);
  }
  
  printf("============================================================\n");
  printf("NSGS Performance Comparison Test\n");
  printf("============================================================\n");
  printf("Data file: %s\n", data_file);
  printf("Number of trials: %d\n", n_trials);
  printf("\n");
  
  /* Load problem from file */
  FrictionContactProblem* problem = frictionContact_new_from_filename(data_file);
  if (!problem) {
    fprintf(stderr, "Error: Failed to load problem from %s\n", data_file);
    return 1;
  }
  
  int nc = problem->numberOfContacts;
  printf("Problem size: %d contacts (%d variables)\n", nc, nc * 3);
  printf("\n");
  
  /* Allocate solution arrays */
  double* reaction = (double*)calloc(nc * 3, sizeof(double));
  double* velocity = (double*)calloc(nc * 3, sizeof(double));
  
  /* Test with different tolerances */
  double tolerances[] = {1e-6, 1e-8, 1e-10, 1e-12};
  int n_tols = sizeof(tolerances) / sizeof(tolerances[0]);
  
  for (int t = 0; t < n_tols; t++) {
    double tol = tolerances[t];
    
    printf("------------------------------------------------------------\n");
    printf("Tolerance: %.2e\n", tol);
    printf("------------------------------------------------------------\n");
    
    /* Setup solver options */
    SolverOptions* options_orig = solver_options_create(SICONOS_FRICTION_3D_NSGS);
    options_orig->dparam[SICONOS_DPARAM_TOL] = tol;
    options_orig->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    solver_options_update_internal(options_orig, 0, SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID);
    options_orig->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = tol * 100;
    options_orig->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    
    SolverOptions* options_gen = solver_options_copy(options_orig);
    
    /* Warmup runs */
    printf("Warmup...\n");
    for (int i = 0; i < 5; i++) {
      int info;
      double* r_warm = (double*)calloc(nc * 3, sizeof(double));
      double* v_warm = (double*)calloc(nc * 3, sizeof(double));
      fc3d_nsgs(problem, r_warm, v_warm, &info, options_orig);
      free(r_warm);
      free(v_warm);
    }
    
    /* Benchmark original fc3d_nsgs */
    SolverBenchmark bench_orig = benchmark_solver(problem, reaction, velocity, 
                                                  options_orig, n_trials, 1, 
                                                  "fc3d_nsgs (original)");
    
    /* Benchmark generic fc3d_nsgs_generic */
    SolverBenchmark bench_gen = benchmark_solver(problem, reaction, velocity,
                                                 options_gen, n_trials, 0,
                                                 "fc3d_nsgs_generic");
    
    /* Print results */
    print_benchmark_results("Original fc3d_nsgs", &bench_orig, n_trials);
    print_benchmark_results("Generic fc3d_nsgs_generic", &bench_gen, n_trials);
    
    /* Compute speedup/overhead */
    double speedup_cpu = bench_orig.cpu_time_avg / bench_gen.cpu_time_avg;
    double speedup_wall = bench_orig.wall_time_avg / bench_gen.wall_time_avg;
    double overhead_cpu = ((bench_gen.cpu_time_avg - bench_orig.cpu_time_avg) / bench_orig.cpu_time_avg) * 100.0;
    
    printf("\nComparison:\n");
    printf("  CPU speedup:   %.2fx (%+.1f%%)\n", speedup_cpu, overhead_cpu);
    printf("  Wall speedup:  %.2fx\n", speedup_wall);
    printf("  Iteration diff: %d\n", bench_gen.iter_avg - bench_orig.iter_avg);
    
    solver_options_delete(options_orig);
    solver_options_delete(options_gen);
  }
  
  /* Final verification with detailed output */
  printf("\n============================================================\n");
  printf("Detailed Verification (tol=1e-8)\n");
  printf("============================================================\n");
  
  SolverOptions* opts = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  opts->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  opts->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  solver_options_update_internal(opts, 0, SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID);
  opts->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  
  double* r_orig = (double*)calloc(nc * 3, sizeof(double));
  double* v_orig = (double*)calloc(nc * 3, sizeof(double));
  double* r_gen = (double*)calloc(nc * 3, sizeof(double));
  double* v_gen = (double*)calloc(nc * 3, sizeof(double));
  
  SolverStats stats_orig = run_solver_once(problem, r_orig, v_orig, opts, 1);
  SolverStats stats_gen = run_solver_once(problem, r_gen, v_gen, opts, 0);
  
  /* Compute solution difference */
  double diff_r = 0.0, diff_v = 0.0;
  for (int i = 0; i < nc * 3; i++) {
    diff_r += (r_orig[i] - r_gen[i]) * (r_orig[i] - r_gen[i]);
    diff_v += (v_orig[i] - v_gen[i]) * (v_orig[i] - v_gen[i]);
  }
  diff_r = sqrt(diff_r);
  diff_v = sqrt(diff_v);
  
  printf("Original:  iter=%d, error=%.6e, time=%.6f ms\n",
         stats_orig.iterations, stats_orig.final_error, stats_orig.cpu_time * 1000);
  printf("Generic:   iter=%d, error=%.6e, time=%.6f ms\n",
         stats_gen.iterations, stats_gen.final_error, stats_gen.cpu_time * 1000);
  printf("Solution diff (reaction): %.6e\n", diff_r);
  printf("Solution diff (velocity): %.6e\n", diff_v);
  
  if (diff_r < 1e-6 && diff_v < 1e-6) {
    printf("\n✓ PASS: Solutions are equivalent\n");
  } else {
    printf("\n✗ FAIL: Solutions differ significantly\n");
  }
  
  printf("============================================================\n");
  
  /* Cleanup */
  free(r_orig);
  free(v_orig);
  free(r_gen);
  free(v_gen);
  free(reaction);
  free(velocity);
  solver_options_delete(opts);
  frictionContactProblem_free(problem);
  
  return 0;
}
