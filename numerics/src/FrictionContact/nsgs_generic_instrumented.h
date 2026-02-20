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

/*!\file nsgs_generic_instrumented.h
 * \brief Instrumented version of nsgs_solve for performance profiling
 */

#ifndef NSGS_GENERIC_INSTRUMENTED_H
#define NSGS_GENERIC_INSTRUMENTED_H

#include "nsgs_generic.h"
#include <sys/time.h>

/* Timer structure for profiling */
typedef struct {
  double total_time;
  double update_local_problem_time;
  double solve_local_time;
  double relaxation_time;
  double accept_local_time;
  double error_computation_time;
  double freezing_time;
  double other_time;
  int update_local_problem_calls;
  int solve_local_calls;
  int relaxation_calls;
  int accept_local_calls;
  int error_computation_calls;
  int freezing_calls;
} NSGSInstrumentedTimers;

/* Global timer instance */
static NSGSInstrumentedTimers nsgs_timers = {0};

/* Reset timers */
static inline void nsgs_timers_reset(void) {
  memset(&nsgs_timers, 0, sizeof(NSGSInstrumentedTimers));
}

/* Print timer results */
static inline void nsgs_timers_print(const char* solver_name, int nb_blocks, int iters) {
  printf("\n============================================================\n");
  printf("Profiling Results: %s\n", solver_name);
  printf("============================================================\n");
  printf("Blocks: %d, Iterations: %d\n", nb_blocks, iters);
  printf("------------------------------------------------------------\n");
  printf("%-30s %12s %10s %10s\n", "Function", "Time (ms)", "Calls", "ms/call");
  printf("------------------------------------------------------------\n");
  
  double total = nsgs_timers.total_time;
  
  #define PRINT_TIMER(name, field) \
    if (nsgs_timers.field > 0) { \
      printf("%-30s %12.4f %10d %10.6f\n", \
             name, \
             nsgs_timers.field * 1000, \
             nsgs_timers.field##_calls, \
             nsgs_timers.field * 1000 / (nsgs_timers.field##_calls > 0 ? nsgs_timers.field##_calls : 1)); \
    }
  
  PRINT_TIMER("update_local_problem", update_local_problem_time);
  PRINT_TIMER("solve_local", solve_local_time);
  PRINT_TIMER("relaxation", relaxation_time);
  PRINT_TIMER("accept_local", accept_local_time);
  PRINT_TIMER("error_computation", error_computation_time);
  PRINT_TIMER("freezing", freezing_time);
  PRINT_TIMER("other", other_time);
  
  #undef PRINT_TIMER
  
  printf("------------------------------------------------------------\n");
  printf("%-30s %12.4f %10s %10s\n", "TOTAL", total * 1000, "-", "-");
  printf("============================================================\n");
  
  /* Print percentages */
  printf("\nTime Breakdown:\n");
  if (nsgs_timers.solve_local_time > 0)
    printf("  solve_local:        %5.1f%%\n", nsgs_timers.solve_local_time / total * 100);
  if (nsgs_timers.update_local_problem_time > 0)
    printf("  update_local:       %5.1f%%\n", nsgs_timers.update_local_problem_time / total * 100);
  if (nsgs_timers.error_computation_time > 0)
    printf("  error_computation:  %5.1f%%\n", nsgs_timers.error_computation_time / total * 100);
  if (nsgs_timers.freezing_time > 0)
    printf("  freezing:           %5.1f%%\n", nsgs_timers.freezing_time / total * 100);
  if (nsgs_timers.relaxation_time > 0)
    printf("  relaxation:         %5.1f%%\n", nsgs_timers.relaxation_time / total * 100);
  if (nsgs_timers.accept_local_time > 0)
    printf("  accept_local:       %5.1f%%\n", nsgs_timers.accept_local_time / total * 100);
  if (nsgs_timers.other_time > 0)
    printf("  other:              %5.1f%%\n", nsgs_timers.other_time / total * 100);
}

/* Get current time in seconds */
static inline double nsgs_get_time(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/** Instrumented NSGS solver - SINGLE LOOP VERSION */
static inline void nsgs_solve_instrumented(void* problem, double* var_z, double* var_x, int* info,
                              SolverOptions* options, NSGSLocalToolkit* toolkit,
                              NSGSProblemData* problem_data) {
  /* Reset timers at start */
  nsgs_timers_reset();
  double global_start = nsgs_get_time();
  
  /* Get solver parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  unsigned int nb_blocks = problem_data->nb_blocks;
  int dim = toolkit->dimension;

  double norm_q = cblas_dnrm2(nb_blocks * dim, problem_data->q, 1);
  double norm_z = 1e24;
  double prev_norm_z = 1e24;

  /* Iteration variables */
  int iter = 0;
  double error = 1.0;
  int hasNotConverged = 1;

  /* Allocate working buffers */
  double* var_z_local = (double*)malloc(dim * sizeof(double));
  double* block_errors = NULL;
  if (toolkit->use_incremental_error || toolkit->use_freezing) {
    block_errors = (double*)malloc(nb_blocks * sizeof(double));
  }

  /* Get local problem */
  void* local_problem = NULL;
  if (toolkit->alloc_local) {
    local_problem = toolkit->alloc_local(problem);
    if (!local_problem) {
      numerics_error("nsgs_solve_instrumented", "Failed to allocate local problem");
      *info = 1;
      goto nsgs_cleanup;
    }
  } else if (toolkit->localproblem) {
    local_problem = toolkit->localproblem;
  }

  /* Allocate block ordering arrays */
  unsigned int* sblocks = NULL;
  unsigned int* freeze_blocks = NULL;

  if (toolkit->use_shuffling) {
    if (toolkit->alloc_shuffled) {
      sblocks = toolkit->alloc_shuffled(problem, options);
    } else {
      sblocks = (unsigned int*)malloc(nb_blocks * sizeof(unsigned int));
      for (unsigned int i = 0; i < nb_blocks; ++i) {
        sblocks[i] = i;
      }
    }
  }

  if (toolkit->use_freezing) {
    if (toolkit->alloc_freezing) {
      freeze_blocks = toolkit->alloc_freezing(problem, options);
    } else {
      freeze_blocks = (unsigned int*)calloc(nb_blocks, sizeof(unsigned int));
    }
  }

  SolverOptions* localsolver_options = NULL;
  if (options->numberOfInternalSolvers > 0) {
    localsolver_options = options->internalSolvers[0];
  }

  if (*info == 0) {
    goto nsgs_cleanup;
  }

  /*****  Main NSGS Iterations *****/
  while ((iter < itermax) && hasNotConverged) {
    ++iter;
    double incremental_error_sum = 0.0;

    /* Set tolerance for internal solver */
    if (toolkit->set_local_tol && localsolver_options) {
      toolkit->set_local_tol(problem, options, localsolver_options, error);
    }

    nsgs_check_freezing_reset(freeze_blocks, nb_blocks, toolkit->use_freezing);
    nsgs_shuffle_blocks(sblocks, nb_blocks, toolkit->use_shuffling ? 2 : 0, iter);

    double tmp_criteria1 = tolerance * tolerance * 100.0 * 100.0;
    double tmp_criteria2 = (prev_norm_z > 0.0) ? 
        (prev_norm_z * prev_norm_z / (nb_blocks * nb_blocks * 1000.0)) : 0.0;

    /* Main loop over blocks */
    for (unsigned int i = 0; i < nb_blocks; ++i) {
      unsigned int block;

      if (sblocks && toolkit->use_shuffling) {
        block = sblocks[i];
      } else {
        block = i;
      }

      if (freeze_blocks && toolkit->use_freezing && freeze_blocks[block] > 0) {
        freeze_blocks[block]--;
        continue;
      }

      /* Copy current solution to local buffer */
      if (dim <= 5) {
        for (int d = 0; d < dim; ++d) {
          var_z_local[d] = var_z[block * dim + d];
        }
      } else {
        cblas_dcopy(dim, &var_z[block * dim], 1, var_z_local, 1);
      }

      /* Update local problem */
      if (toolkit->update_local_problem) {
        double t_start = nsgs_get_time();
        toolkit->update_local_problem(block, problem, local_problem, var_z, options);
        nsgs_timers.update_local_problem_time += nsgs_get_time() - t_start;
        nsgs_timers.update_local_problem_calls++;
      }

      /* Solve local problem */
      if (toolkit->solve_local) {
        double t_start = nsgs_get_time();
        toolkit->solve_local(local_problem, localsolver_options, var_z_local, NULL);
        nsgs_timers.solve_local_time += nsgs_get_time() - t_start;
        nsgs_timers.solve_local_calls++;
      }

      /* Perform relaxation */
      if (toolkit->use_relaxation && toolkit->relaxation) {
        double t_start = nsgs_get_time();
        toolkit->relaxation(var_z_local, &var_z[block * dim], toolkit->omega);
        nsgs_timers.relaxation_time += nsgs_get_time() - t_start;
        nsgs_timers.relaxation_calls++;
      }

      /* Compute incremental error and freezing criteria */
      double t_start_error = nsgs_get_time();
      double incr_err_2 = 0.0;
      double sq_norm_local = 0.0;
      if (toolkit->use_incremental_error || toolkit->use_freezing) {
        for (int d = 0; d < dim; ++d) {
          double diff = var_z_local[d] - var_z[block * dim + d];
          incr_err_2 += diff * diff;
          sq_norm_local += var_z_local[d] * var_z_local[d];
        }
        if (toolkit->use_incremental_error) {
          incremental_error_sum += incr_err_2;
          block_errors[block] = incr_err_2;
        }
      }
      nsgs_timers.error_computation_time += nsgs_get_time() - t_start_error;
      nsgs_timers.error_computation_calls++;

      /* Check freezing */
      double t_start_freeze = nsgs_get_time();
      if (freeze_blocks && toolkit->use_freezing && toolkit->freezing_iter > 0 && iter >= 10) {
        int relative_conv = incr_err_2 <= tmp_criteria1 * sq_norm_local;
        int small_sol = sq_norm_local <= tmp_criteria2;
        if (relative_conv || small_sol) {
          freeze_blocks[block] = toolkit->freezing_iter;
        }
      }
      nsgs_timers.freezing_time += nsgs_get_time() - t_start_freeze;
      nsgs_timers.freezing_calls++;

      /* Accept local solution */
      double t_start_accept = nsgs_get_time();
      if (toolkit->accept_local) {
        toolkit->accept_local(local_problem, localsolver_options, block, iter, var_z,
                              var_z_local);
      } else {
        if (dim <= 5) {
          for (int d = 0; d < dim; ++d) {
            var_z[block * dim + d] = var_z_local[d];
          }
        } else {
          cblas_dcopy(dim, var_z_local, 1, &var_z[block * dim], 1);
        }
      }
      nsgs_timers.accept_local_time += nsgs_get_time() - t_start_accept;
      nsgs_timers.accept_local_calls++;
    }

    /* Compute global error */
    prev_norm_z = norm_z;
    norm_z = cblas_dnrm2(nb_blocks * dim, var_z, 1);
    
    double incremental_error = 0.0;
    if (toolkit->use_incremental_error) {
      incremental_error = sqrt(incremental_error_sum);
      if (norm_z > 0.0) {
        incremental_error /= norm_z;
      }
    }

    hasNotConverged = 1;
    error = incremental_error;
    double full_error = 0.0;
    int computed_full_error = 0;

    nsgs_early_tolerance_adapt(problem, var_z, var_x, options, toolkit,
                                incremental_error, &tolerance, &full_error, iter);
    if (iter == 10) computed_full_error = 1;

    if (incremental_error <= tolerance) {
      hasNotConverged = nsgs_check_convergence_with_full_error(
          problem, var_z, var_x, options, toolkit, incremental_error, &tolerance,
          &full_error, computed_full_error, localsolver_options, iter);
      error = full_error;
    } else {
      if (toolkit->check_convergence) {
        hasNotConverged = toolkit->check_convergence(incremental_error, tolerance, iter, options);
      } else {
        hasNotConverged = nsgs_determine_convergence(incremental_error, tolerance, iter, options);
      }
    }

    nsgs_print_iteration_stats(iter, incremental_error, full_error, tolerance,
                                toolkit->user_tolerance, 
                                nsgs_count_frozen_percent(freeze_blocks, nb_blocks, toolkit->use_freezing),
                                hasNotConverged, toolkit->verbose);

    if (toolkit->stats_callback) {
      toolkit->stats_callback(problem, options, var_z, var_x, error);
    }
  }

  if (toolkit->verbose > 0 && iter > 0) {
    printf("\n");
  }

  if (toolkit->compute_error) {
    double final_full_error;
    int final_converged = !nsgs_check_full_error_convergence(
        problem, var_z, var_x, options, toolkit, toolkit->user_tolerance, &final_full_error);
    
    if (toolkit->verbose > 0) {
      printf("Final check: Full error = %.6e, User tolerance = %.6e, Converged = %s\n", 
             final_full_error, toolkit->user_tolerance, final_converged ? "YES" : "NO");
    }
    
    error = final_full_error;
    *info = final_converged ? 0 : 1;
    hasNotConverged = final_converged ? 0 : 1;
  }

  if (error <= tolerance) {
    *info = 0;
  } else {
    *info = 1;
  }

  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  dparam[SICONOS_DPARAM_RESIDU] = error;

nsgs_cleanup:
  free(var_z_local);
  if (block_errors) free(block_errors);
  if (sblocks) free(sblocks);
  if (freeze_blocks) free(freeze_blocks);
  
  /* Record total time */
  nsgs_timers.total_time = nsgs_get_time() - global_start;
  /* Calculate other time */
  nsgs_timers.other_time = nsgs_timers.total_time 
                         - nsgs_timers.update_local_problem_time
                         - nsgs_timers.solve_local_time
                         - nsgs_timers.relaxation_time
                         - nsgs_timers.accept_local_time
                         - nsgs_timers.error_computation_time
                         - nsgs_timers.freezing_time;
}

#endif /* NSGS_GENERIC_INSTRUMENTED_H */
