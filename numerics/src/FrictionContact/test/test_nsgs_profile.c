/* Test to profile nsgs_solve vs fc3d_nsgs with detailed timing */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "NumericsFwd.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#include "numerics_verbose.h"
#include "nsgs_generic_instrumented.h"

/* External timer struct from instrumented header */
extern NSGSInstrumentedTimers nsgs_timers;

/* Forward declaration of instrumented fc3d_nsgs - we'll need to add this to the library */
/* For now, we'll just profile the generic solver */

/* Get wall clock time in seconds */
static double get_wall_time(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Wrapper to call instrumented nsgs_solve from fc3d_nsgs_generic */
static void fc3d_nsgs_generic_instrumented(FrictionContactProblem* problem, 
                                           double* reaction, double* velocity,
                                           int* info, SolverOptions* options) {
  if (!problem || !reaction || !info || !options) {
    numerics_error("fc3d_nsgs_generic_instrumented", "Invalid input arguments");
    return;
  }
  
  if (options->numberOfInternalSolvers < 1) {
    numerics_error("fc3d_nsgs_generic_instrumented",
                   "The NSGS method needs options for the internal solvers");
    *info = 1;
    return;
  }

  /* Extract diagonal blocks */
  void* original_matrix1 = nsgs_generic_extract_diagonal_blocks(problem->M, problem->dimension);

  /* Create local problem */
  SolverPtr local_solver = NULL;
  FrictionContactProblem* localproblem = fc3d_nsgs_local_problem_allocate(problem);
  if (!localproblem) {
    *info = 1;
    return;
  }
  
  /* Initialize local solver based on type */
  SolverOptions* localsolver_options = options->internalSolvers[0];
  switch (localsolver_options->solverId) {
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone:
      fc3d_projection_initialize(problem, localproblem);
      local_solver = &fc3d_projectionOnCone_solve;
      break;
    case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID:
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
      local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      break;
    default:
      fc3d_local_problem_free(localproblem, problem);
      *info = 1;
      return;
  }
  
  localsolver_options->solverData = (void*)local_solver;

  /* Setup problem data */
  NSGSProblemData problem_data = {
    .nb_blocks = problem->numberOfContacts,
    .q = problem->q,
    .M = problem->M,
    .mu = problem->mu,
    .mu_r = NULL,
    .storage_type = problem->M->storageType,
    .dimension = 3
  };

  /* Setup toolkit */
  NSGSLocalToolkit toolkit = {
    .update_local_problem = fc3d_nsgs_update_wrapper,
    .solve_local = fc3d_nsgs_solve_local_wrapper,
    .compute_error = fc3d_nsgs_compute_error_wrapper,
    .incremental_error = fc3d_incr_error_3,
    .accept_local = fc3d_nsgs_accept_local_wrapper,
    .check_convergence = NULL,
    .alloc_local = NULL,
    .set_local_tol = fc3d_nsgs_set_tol_wrapper,
    .stats_callback = NULL,
    .relaxation = fc3d_relaxation_3,
    .squared_norm = fc3d_squared_norm_3,
    .should_freeze = NULL,
    .alloc_freezing = NULL,
    .alloc_shuffled = NULL,
    .localproblem = localproblem,
    .verbose = 0,  /* Disable verbose output during profiling */
    .user_tolerance = options->dparam[SICONOS_DPARAM_TOL],
    .error_eval_freq = 0,
    .dimension = 3,
    .omega = options->dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE],
    .use_freezing = options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0,
    .use_shuffling = options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] != SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE,
    .use_relaxation = options->iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE,
    .use_incremental_error = 1,
    .filter_local_sol = options->iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION],
    .error_eval_type = options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION],
    .freezing_iter = options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]
  };

  /* Call instrumented solver */
  nsgs_solve_instrumented(problem, reaction, velocity, info, options, &toolkit, &problem_data);

  /* Cleanup */
  fc3d_nsgs_local_problem_free(localproblem, problem);
  localsolver_options->solverData = NULL;
  nsgs_generic_free_diagonal_blocks(problem->M, original_matrix1);
}

/* Run and profile a solver */
typedef struct {
  double total_time;
  int iterations;
  int info;
  double final_error;
} ProfileResult;

static ProfileResult profile_solver(FrictionContactProblem* problem,
                                    double* reaction, double* velocity,
                                    SolverOptions* options,
                                    int use_original) {
  ProfileResult result = {0};
  int nc = problem->numberOfContacts;
  
  double* r_tmp = (double*)calloc(nc * 3, sizeof(double));
  double* v_tmp = (double*)calloc(nc * 3, sizeof(double));
  
  double start = get_wall_time();
  
  if (use_original) {
    fc3d_nsgs(problem, r_tmp, v_tmp, &result.info, options);
  } else {
    fc3d_nsgs_generic_instrumented(problem, r_tmp, v_tmp, &result.info, options);
  }
  
  result.total_time = get_wall_time() - start;
  result.iterations = options->iparam[SICONOS_IPARAM_ITER_DONE];
  result.final_error = options->dparam[SICONOS_DPARAM_RESIDU];
  
  /* Copy solution back */
  memcpy(reaction, r_tmp, nc * 3 * sizeof(double));
  memcpy(velocity, v_tmp, nc * 3 * sizeof(double));
  
  free(r_tmp);
  free(v_tmp);
  
  return result;
}

int main(int argc, char** argv) {
  printf("============================================================\n");
  printf("NSGS Profiling Test\n");
  printf("============================================================\n");
  
  const char* data_file = (argc > 1) ? argv[1] : "./data/FC3D_Example1.dat";
  
  FrictionContactProblem* problem = frictionContact_new_from_filename(data_file);
  if (!problem) {
    fprintf(stderr, "Failed to load %s\n", data_file);
    return 1;
  }
  
  int nc = problem->numberOfContacts;
  printf("Problem: %s (%d contacts, %d variables)\n\n", data_file, nc, nc*3);
  
  /* Setup options */
  SolverOptions* opts_orig = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  opts_orig->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  opts_orig->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  solver_options_update_internal(opts_orig, 0, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  opts_orig->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  
  SolverOptions* opts_gen = solver_options_copy(opts_orig);
  
  /* Allocate solution arrays */
  double* r_orig = (double*)calloc(nc * 3, sizeof(double));
  double* v_orig = (double*)calloc(nc * 3, sizeof(double));
  double* r_gen = (double*)calloc(nc * 3, sizeof(double));
  double* v_gen = (double*)calloc(nc * 3, sizeof(double));
  
  /* Profile original fc3d_nsgs */
  printf("Profiling ORIGINAL fc3d_nsgs...\n");
  printf("------------------------------------------------------------\n");
  ProfileResult res_orig = profile_solver(problem, r_orig, v_orig, opts_orig, 1);
  printf("Total time:   %.4f ms\n", res_orig.total_time * 1000);
  printf("Iterations:   %d\n", res_orig.iterations);
  printf("Final error:  %.6e\n", res_orig.final_error);
  printf("Converged:    %s\n\n", res_orig.info == 0 ? "YES" : "NO");
  
  /* Profile generic fc3d_nsgs_generic with instrumentation */
  printf("Profiling GENERIC fc3d_nsgs_generic (instrumented)...\n");
  printf("------------------------------------------------------------\n");
  ProfileResult res_gen = profile_solver(problem, r_gen, v_gen, opts_gen, 0);
  
  /* Print instrumented results */
  nsgs_timers_print("fc3d_nsgs_generic", nc, res_gen.iterations);
  
  printf("\nComparison:\n");
  printf("------------------------------------------------------------\n");
  printf("Total time (original):  %.4f ms\n", res_orig.total_time * 1000);
  printf("Total time (generic):   %.4f ms\n", res_gen.total_time * 1000);
  printf("Overhead:               %.1f%%\n", 
         (res_gen.total_time - res_orig.total_time) / res_orig.total_time * 100);
  printf("Iterations (original):  %d\n", res_orig.iterations);
  printf("Iterations (generic):   %d\n", res_gen.iterations);
  
  /* Compute solution difference */
  double diff = 0.0;
  for (int i = 0; i < nc * 3; i++) {
    double d = r_orig[i] - r_gen[i];
    diff += d * d;
  }
  diff = sqrt(diff);
  printf("Solution difference:    %.6e\n", diff);
  
  /* Cleanup */
  free(r_orig); free(v_orig);
  free(r_gen); free(v_gen);
  solver_options_delete(opts_orig);
  solver_options_delete(opts_gen);
  frictionContactProblem_free(problem);
  
  return 0;
}

/* Wrapper functions needed by the instrumented solver */
extern void fc3d_nsgs_update(int, FrictionContactProblem*, FrictionContactProblem*, 
                             double*, SolverOptions*);
extern double fc3d_compute_error;

static void fc3d_nsgs_update_wrapper(unsigned int block, void* problem, 
                                     void* local_problem, double* var_z, 
                                     SolverOptions* options) {
  fc3d_nsgs_update(block, (FrictionContactProblem*)problem,
                   (FrictionContactProblem*)local_problem, var_z, options);
}

static void fc3d_nsgs_solve_local_wrapper(void* local_problem, 
                                          SolverOptions* localsolver_options,
                                          double* var_z_local, 
                                          double* localsolver_options_data) {
  (void)localsolver_options_data;
  SolverPtr local_solver = (SolverPtr)localsolver_options->solverData;
  if (local_solver && local_problem) {
    local_solver((FrictionContactProblem*)local_problem, var_z_local, localsolver_options);
  }
}

static double fc3d_nsgs_compute_error_wrapper(void* problem, double* var_z,
                                              double* var_x, SolverOptions* options) {
  double error = 0.0;
  double norm_q = 0.0;
  fc3d_compute_error((FrictionContactProblem*)problem, var_z, var_x,
                     options->dparam[SICONOS_DPARAM_TOL], options, norm_q, &error);
  return error;
}

static void fc3d_nsgs_set_tol_wrapper(void* problem, SolverOptions* options,
                                      SolverOptions* localsolver_options, double error) {
  (void)problem;
  (void)options;
  (void)error;
  /* Default: no tolerance setting */
}

static void fc3d_nsgs_accept_local_wrapper(void* local_problem, SolverOptions* options,
                                           unsigned int block, int iter,
                                           double* var_z_global, double* var_z_local) {
  (void)local_problem;
  (void)iter;
  
  double local_residual = options->dparam[SICONOS_DPARAM_RESIDU];
  if (isnan(local_residual) || isinf(local_residual) || local_residual > 1.0) {
    numerics_printf("Discard local reaction for block %i at iteration %i with local_error = %e",
                    block, iter, local_residual);
    return;
  }
  
  var_z_global[block * 3 + 0] = var_z_local[0];
  var_z_global[block * 3 + 1] = var_z_local[1];
  var_z_global[block * 3 + 2] = var_z_local[2];
}

static void fc3d_relaxation_3(double* z_new, double* z_old, double omega) {
  z_new[0] = omega * z_new[0] + (1.0 - omega) * z_old[0];
  z_new[1] = omega * z_new[1] + (1.0 - omega) * z_old[1];
  z_new[2] = omega * z_new[2] + (1.0 - omega) * z_old[2];
}

static double fc3d_squared_norm_3(double* z) {
  return z[0]*z[0] + z[1]*z[1] + z[2]*z[2];
}

static double fc3d_incr_error_3(double* z_new, double* z_old) {
  double diff[3] = {z_new[0] - z_old[0], z_new[1] - z_old[1], z_new[2] - z_old[2]};
  return diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
}
