/* Detailed comparison test between fc3d_nsgs and fc3d_nsgs_generic */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "NumericsFwd.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
/* #include "numerics_verbose.h" */

/* Get CPU time in seconds */
static double get_cpu_time(void) {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec * 1e-6 +
         (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec * 1e-6;
}

/* Get wall clock time */
static double get_wall_time(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Result structure */
typedef struct {
  double cpu_time;
  double wall_time;
  int iterations;
  int info;
  double final_error;
  double reaction_norm;
  double velocity_norm;
} SolverResult;

/* Run solver and collect results */
static SolverResult run_solver(FrictionContactProblem* problem,
                                double* reaction, double* velocity,
                                SolverOptions* options, int use_original) {
  SolverResult res = {0};
  int nc = problem->numberOfContacts;
  
  /* Reset initial guess */
  memset(reaction, 0, nc * 3 * sizeof(double));
  memset(velocity, 0, nc * 3 * sizeof(double));
  
  /* IMPORTANT: Reset iteration counter before each run */
  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
  
  /* Initialize info to non-zero to force solver to run (not early exit) */
  res.info = -1;
  
  double cpu_start = get_cpu_time();
  double wall_start = get_wall_time();
  
  if (use_original) {
    fc3d_nsgs(problem, reaction, velocity, &res.info, options);
  } else {
    fc3d_nsgs_generic(problem, reaction, velocity, &res.info, options);
  }
  
  res.cpu_time = get_cpu_time() - cpu_start;
  res.wall_time = get_wall_time() - wall_start;
  res.iterations = options->iparam[SICONOS_IPARAM_ITER_DONE];
  res.final_error = options->dparam[SICONOS_DPARAM_RESIDU];
  
  /* Compute norms */
  for (int i = 0; i < nc * 3; i++) {
    res.reaction_norm += reaction[i] * reaction[i];
    res.velocity_norm += velocity[i] * velocity[i];
  }
  res.reaction_norm = sqrt(res.reaction_norm);
  res.velocity_norm = sqrt(res.velocity_norm);
  
  return res;
}

/* Compare two solutions */
typedef struct {
  double reaction_diff;
  double velocity_diff;
  double max_reaction_diff;
  int max_diff_index;
} SolutionDiff;

static SolutionDiff compare_solutions(int n, double* r1, double* v1, double* r2, double* v2) {
  SolutionDiff diff = {0};
  diff.max_reaction_diff = 0.0;
  diff.max_diff_index = -1;
  
  for (int i = 0; i < n; i++) {
    double dr = r1[i] - r2[i];
    double dv = v1[i] - v2[i];
    diff.reaction_diff += dr * dr;
    diff.velocity_diff += dv * dv;
    
    if (fabs(dr) > diff.max_reaction_diff) {
      diff.max_reaction_diff = fabs(dr);
      diff.max_diff_index = i;
    }
  }
  diff.reaction_diff = sqrt(diff.reaction_diff);
  diff.velocity_diff = sqrt(diff.velocity_diff);
  
  return diff;
}

/* Print comparison table */
static void print_comparison(const char* label, 
                              SolverResult* orig, SolverResult* gen,
                              SolutionDiff* diff) {
  printf("\n╔══════════════════════════════════════════════════════════════════╗\n");
  printf("║  %s\n", label);
  printf("╠══════════════════════════════════════════════════════════════════╣\n");
  printf("║  Metric                │  Original      │  Generic       │ Diff   ║\n");
  printf("╠══════════════════════════════════════════════════════════════════╣\n");
  printf("║  CPU Time (ms)         │  %12.4f  │  %12.4f  │ %+6.1f%% ║\n",
         orig->cpu_time * 1000, gen->cpu_time * 1000,
         (gen->cpu_time - orig->cpu_time) / orig->cpu_time * 100);
  printf("║  Wall Time (ms)        │  %12.4f  │  %12.4f  │ %+6.1f%% ║\n",
         orig->wall_time * 1000, gen->wall_time * 1000,
         (gen->wall_time - orig->wall_time) / orig->wall_time * 100);
  printf("║  Iterations            │  %12d  │  %12d  │ %+6d  ║\n",
         orig->iterations, gen->iterations, gen->iterations - orig->iterations);
  printf("║  Final Error           │  %12.4e  │  %12.4e  │        ║\n",
         orig->final_error, gen->final_error);
  printf("║  Converged             │  %12s  │  %12s  │        ║\n",
         orig->info == 0 ? "YES" : "NO", gen->info == 0 ? "YES" : "NO");
  printf("║  Reaction Norm         │  %12.4e  │  %12.4e  │        ║\n",
         orig->reaction_norm, gen->reaction_norm);
  printf("╠══════════════════════════════════════════════════════════════════╣\n");
  printf("║  Solution Differences                                            ║\n");
  printf("║  Reaction L2-norm:     %12.4e                                 ║\n", diff->reaction_diff);
  printf("║  Velocity L2-norm:     %12.4e                                 ║\n", diff->velocity_diff);
  printf("║  Max reaction diff:    %12.4e  at index %d                     ║\n",
         diff->max_reaction_diff, diff->max_diff_index);
  printf("╚══════════════════════════════════════════════════════════════════╝\n");
}

/* Test on a single problem */
static void test_problem(const char* filename, double tol) {
  printf("\n");
  printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  printf("  Problem: %s (tol=%.2e)\n", filename, tol);
  printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  
  FrictionContactProblem* problem = frictionContact_new_from_filename(filename);
  if (!problem) {
    printf("  ERROR: Failed to load %s\n", filename);
    return;
  }
  
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  printf("  Size: %d contacts (%d variables)\n\n", nc, n);
  
  /* Setup solver options - separate for each solver */
  SolverOptions* opts_orig = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  opts_orig->dparam[SICONOS_DPARAM_TOL] = tol;
  opts_orig->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  solver_options_update_internal(opts_orig, 0, SICONOS_ONECONE_NSN_GP_HYBRID);
  opts_orig->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = tol * 100;
  opts_orig->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  
  /* Copy for generic solver */
  SolverOptions* opts_gen = solver_options_copy(opts_orig);
  
  /* Allocate solution arrays */
  double *r_orig = (double*)calloc(n, sizeof(double));
  double *v_orig = (double*)calloc(n, sizeof(double));
  double *r_gen = (double*)calloc(n, sizeof(double));
  double *v_gen = (double*)calloc(n, sizeof(double));
  
  /* Run original solver */
  printf("  Running original fc3d_nsgs...\n");
  SolverResult res_orig = run_solver(problem, r_orig, v_orig, opts_orig, 1);
  
  /* Run generic solver */
  printf("  Running generic fc3d_nsgs_generic...\n");
  SolverResult res_gen = run_solver(problem, r_gen, v_gen, opts_gen, 0);
  
  /* Compare solutions */
  SolutionDiff diff = compare_solutions(n, r_orig, v_orig, r_gen, v_gen);
  
  /* Print results */
  print_comparison(filename, &res_orig, &res_gen, &diff);
  
  /* First 6 reaction values */
  printf("\n  First 6 reaction values:\n");
  printf("  %-6s  %-16s  %-16s  %-12s\n", "Index", "Original", "Generic", "Diff");
  printf("  %-6s  %-16s  %-16s  %-12s\n", "------", "----------------", "----------------", "------------");
  for (int i = 0; i < 6 && i < n; i++) {
    printf("  %-6d  %16.6e  %16.6e  %12.4e\n", 
           i, r_orig[i], r_gen[i], fabs(r_orig[i] - r_gen[i]));
  }
  
  /* Cleanup */
  free(r_orig); free(v_orig);
  free(r_gen); free(v_gen);
  solver_options_delete(opts_orig);
  solver_options_delete(opts_gen);
  frictionContactProblem_free(problem);
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  
  printf("\n");
  printf("╔══════════════════════════════════════════════════════════════════╗\n");
  printf("║         NSGS Detailed Comparison Test                            ║\n");
  printf("║         Original fc3d_nsgs vs Generic fc3d_nsgs_generic          ║\n");
  printf("╚══════════════════════════════════════════════════════════════════╝\n");
  
  /* Test with different tolerances on different problems */
  struct {
    const char* file;
    double tol;
  } tests[] = {
    {"./data/FC3D_Example1.dat", 1e-6},
    {"./data/FC3D_Example1.dat", 1e-8},
    {"./data/Confeti-ex13-Fc3D-SBM.dat", 1e-6},
    {"./data/Confeti-ex13-Fc3D-SBM.dat", 1e-8},
    {"./data/KaplasTower-i1061-4.hdf5.dat", 1e-6},
  };
  
  int n_tests = sizeof(tests) / sizeof(tests[0]);
  
  for (int i = 0; i < n_tests; i++) {
    test_problem(tests[i].file, tests[i].tol);
  }
  
  printf("\n");
  printf("╔══════════════════════════════════════════════════════════════════╗\n");
  printf("║         Summary                                                  ║\n");
  printf("╚══════════════════════════════════════════════════════════════════╝\n");
  printf("\nKey Observations:\n");
  printf("  1. Iteration counts should be similar or identical\n");
  printf("  2. CPU time shows overhead of generic implementation\n");
  printf("  3. Solution differences indicate numerical equivalence\n");
  printf("  4. Convergence behavior should match closely\n");
  printf("\n");
  
  return 0;
}
