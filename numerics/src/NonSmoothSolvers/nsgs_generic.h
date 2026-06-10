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

/*!\file nsgs_generic.h
 * \brief Generic Nonsmooth Gauss-Seidel (NSGS) solver framework
 *
 * This header provides a generic implementation of the NSGS algorithm
 * that can be used for various friction contact problems (2D, 3D, rolling, etc.)
 * by providing appropriate callback functions.
 *
 * Terminology:
 * - "block" : a subsystem (e.g., a contact in friction problems)
 * - "var_z" : primal variable (reaction, lambda, etc.)
 * - "var_x" : dual variable (velocity, relative velocity, etc.)
 */

#ifndef NSGS_GENERIC_H
#define NSGS_GENERIC_H

#include <float.h>  // for DBL_EPSILON
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NumericsMatrix.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "SparseBlockMatrix.h" /* for SparseBlockStructuredMatrix */
#include "numerics_verbose.h"

/** Callback type for updating local problem for a given block */
typedef int (*NSGSUpdateLocalProblem)(unsigned int block, void* problem, void* local_problem,
                                       double* var_z_global, SolverOptions* options);

/** Callbavk type for solving the local problem */
typedef int (*NSGSSolveLocal)(void* local_problem, SolverOptions* options,
                               double* var_z_local, double* localsolver_options_data);

/** Callback type for computing global error */
typedef double (*NSGSComputeError)(void* problem, double* var_z, double* var_x,
                                   SolverOptions* options);

/** Callback type for copying old local solution from global to local buffer (warm start) */
typedef void (*NSGSCopyLocal)(unsigned int block, double* var_z_global, double* var_z_local);

/** Callback type for computing incremental error (squared) between new and old local solution
 */
typedef double (*NSGSIncrError)(double* var_z_local_new, double* var_z_local_old);

/** Callback type for accepting local solution into global solution */
typedef int (*NSGSAcceptLocal)(void* local_problem, SolverOptions* options,
                                unsigned int block, int iter, double* var_z_global,
                                double* var_z_local);

/** Callback type for checking convergence */
typedef int (*NSGSCheckConvergence)(double error, double tolerance, int iter,
                                    SolverOptions* options);

/** Callback type for allocating/freezing blocks array */
typedef unsigned int* (*NSGSAllocFreezing)(void* problem, SolverOptions* options);

/** Callback type for allocating shuffled blocks array */
typedef unsigned int* (*NSGSAllocShuffled)(void* problem, SolverOptions* options);

/** Callback type for allocating local problem */
typedef void* (*NSGSAllocLocalProblem)(void* problem);

/** Callback type for setting internal solver tolerance */
typedef void (*NSGSSetLocalTol)(void* problem, SolverOptions* options,
                                SolverOptions* localsolver_options, double error);

/** Callback type for statistics callback */
typedef void (*NSGSStatsCallback)(void* problem, SolverOptions* options, double* var_z,
                                  double* var_x, double error);

/** Callback type for performing relaxation on local solution */
typedef void (*NSGSRelaxation)(double* var_z_local_new, double* var_z_local_old, double omega);

/** Callback type for computing squared norm of a local solution */
typedef double (*NSGSSquaredNorm)(double* var_z_local);

/** Callback type for checking if a block should be frozen */
typedef int (*NSGSShouldFreeze)(double incremental_error, double squared_norm_local,
                                double tolerance, double norm_z_global, unsigned int nb_blocks,
                                int iter);

/** Generic local problem toolkit structure
 *
 * This structure holds all the callbacks needed to specialize the NSGS
 * algorithm for a specific problem type.
 */
typedef struct {
  NSGSUpdateLocalProblem update_local_problem;
  NSGSCopyLocal copy_local;
  NSGSSolveLocal solve_local;
  NSGSComputeError compute_error;
  NSGSIncrError incremental_error;
  NSGSAcceptLocal accept_local;
  NSGSCheckConvergence check_convergence;
  NSGSAllocFreezing alloc_freezing;
  NSGSAllocShuffled alloc_shuffled;
  NSGSAllocLocalProblem alloc_local;
  NSGSSetLocalTol set_local_tol;
  NSGSStatsCallback stats_callback;
  NSGSRelaxation relaxation;
  NSGSSquaredNorm squared_norm;
  NSGSShouldFreeze should_freeze;

  /* Problem-specific data */
  int dimension;             /* Dimension per block (2, 3, 5, etc.) */
  double omega;              /* Relaxation parameter */
  int use_freezing;          /* Whether to use block freezing */
  int use_shuffling;         /* Whether to shuffle blocks */
  int use_incremental_error; /* Whether to use incremental error for convergence */
  int use_relaxation;        /* Whether to use relaxation */
  int filter_local_sol;      /* Whether to filter local solutions */
  int error_eval_type;   /* Error evaluation mode (0=FULL, 1=LIGHT_WITH_FULL_FINAL, 2=LIGHT,
                            3=ADAPTIVE) */
  int error_eval_freq;   /* Full error evaluation frequency (0=every iter, N=every N iter) */
  int freezing_iter;     /* Number of iterations to freeze a block */
  int verbose;           /* Verbosity level (0=none, 1=iterations, 2=detailed) */
  double user_tolerance; /* User-requested tolerance (for reference, adaptive logic uses
                            working tolerance) */
  void* localproblem;    /* Pointer to local problem (passed to callbacks) */
} NSGSLocalToolkit;

/** Generic problem data structure
 *
 * This holds the common data needed by the NSGS algorithm.
 */
typedef struct {
  unsigned int nb_blocks; /* Number of blocks */
  double* q;              /* Right-hand side vector */
  void* M;                /* Matrix (NumericsMatrix*) */
  double* mu;             /* Friction coefficients (optional) */
  double* mu_r;           /* Rolling friction coefficients (optional) */
  int storage_type;       /* Matrix storage type */
  int dimension;          /* Problem dimension per block */
} NSGSProblemData;

/** Generic local problem structure for block-structured matrices
 *
 * This structure holds the local problem data for a single block.
 * It assumes all blocks have the same dimension.
 */
typedef struct {
  double* M_local;    /* Local matrix (dimension x dimension) */
  double* q_local;    /* Local RHS vector (dimension) */
  double* mu_local;   /* Local friction coefficient (optional) */
  double* mu_r_local; /* Local rolling friction (optional) */
  int dimension;      /* Dimension of this local problem */
} NSGSLocalProblem;

/** Generic local problem update function
 *
 * Updates the local problem for a given block by:
 * 1. Extracting the diagonal block from the global matrix M
 * 2. Computing local q = q_global[block] + off-diagonal contributions
 *
 * This function works with any block-structured matrix where all blocks
 * have the same dimension. For common dimensions (2, 3, 5), it uses optimized
 * BLAS-based functions. For other dimensions, it falls back to generic loops.
 *
 * \param[in] block Block index
 * \param[in] problem_data Problem data containing M, q, etc.
 * \param[in,out] local_problem Local problem structure to update
 * \param[in] var_z_global Global solution vector
 * \param[in] dim Dimension per block
 */
static inline void nsgs_generic_update_local_problem(unsigned int block,
                                                     NSGSProblemData* problem_data,
                                                     NSGSLocalProblem* local_problem,
                                                     double* var_z_global, int dim) {
  NumericsMatrix* M = (NumericsMatrix*)problem_data->M;
  double* q = problem_data->q;
  int block_start = block * dim;

  /* =====================================================================
   * Step 1: Extract diagonal block into local M
   * ===================================================================== */

  /* Check for pre-extracted diagonal blocks (NM_SPARSE with matrix1) */
  if (M->storageType == NM_SPARSE && M->matrix1) {
    /* Use pre-extracted block from SparseBlockStructuredMatrix */
    local_problem->M_local = M->matrix1->block[block];
  } else {
    /* Extract diagonal block on-the-fly */
    switch (dim) {
      case 2:
        NM_extract_diag_block2(M, block, &local_problem->M_local);
        break;
      case 3:
        NM_extract_diag_block3(M, block, &local_problem->M_local);
        break;
      case 5:
        NM_extract_diag_block5(M, block, &local_problem->M_local);
        break;
      default:
        /* Generic fallback for other dimensions */
        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            local_problem->M_local[i * dim + j] =
                NM_get_value(M, block_start + i, block_start + j);
          }
        }
        break;
    }
  }

  /* =====================================================================
   * Step 2: Compute local q = q_block + off-diagonal contributions
   * ===================================================================== */
  /* First copy the local q from global q */
  for (int i = 0; i < dim; ++i) {
    local_problem->q_local[i] = q[block_start + i];
  }

  /* Add off-diagonal contributions using dimension-specific optimized functions */
  switch (dim) {
    case 2:
      NM_row_prod_no_diag2(problem_data->nb_blocks * 2, block, block_start, M, var_z_global,
                           local_problem->q_local, false);
      break;
    case 3:
      NM_row_prod_no_diag3(problem_data->nb_blocks * 3, block, block_start, M, var_z_global,
                           local_problem->q_local, false);
      break;
    default:
      /* Generic fallback for other dimensions using NM_row_prod_no_diag */
      NM_row_prod_no_diag(problem_data->nb_blocks * dim, dim, block, block_start, M,
                          var_z_global, local_problem->q_local, NULL, false);
      break;
  }

  /* Copy friction coefficients if present */
  if (problem_data->mu && local_problem->mu_local) {
    local_problem->mu_local[0] = problem_data->mu[block];
  }
  if (problem_data->mu_r && local_problem->mu_r_local) {
    local_problem->mu_r_local[0] = problem_data->mu_r[block];
  }
}

/** Wrapper for generic local problem update with callback signature
 *
 * This wrapper adapts nsgs_generic_update_local_problem to the NSGSUpdateLocalProblem
 * callback signature, storing the NSGSProblemData pointer in options->solverData.
 *
 * Usage:
 *   - Set options->solverData = &problem_data before calling nsgs_solve
 *   - Set toolkit->update_local_problem = nsgs_generic_update_local_problem_callback
 *   - The local_problem pointer should point to an NSGSLocalProblem structure
 *
 * \param[in] block Block index
 * \param[in] problem Problem pointer (unused, data comes from options->solverData)
 * \param[in,out] local_problem Local problem structure (NSGSLocalProblem*)
 * \param[in] var_z_global Global solution vector
 * \param[in] options Solver options (contains problem_data in solverData)
 */
static inline int nsgs_generic_update_local_problem_callback(unsigned int block,
                                                              void* problem,
                                                              void* local_problem,
                                                              double* var_z_global,
                                                              SolverOptions* options) {
  (void)problem; /* Not used, data comes from options->solverData */

  NSGSProblemData* problem_data = (NSGSProblemData*)options->solverData;
  if (!problem_data) {
    return numerics_error("nsgs_generic_update_local_problem_callback",
                   "options->solverData must contain NSGSProblemData pointer");
  }

  nsgs_generic_update_local_problem(block, problem_data, (NSGSLocalProblem*)local_problem,
                                    var_z_global, problem_data->dimension);
  return 0;
}

/** Extract diagonal blocks from matrix for efficient access
 *
 * For NM_SPARSE matrices, extracts diagonal blocks and stores them in M->matrix1
 * as a SparseBlockStructuredMatrix. This allows O(1) access to diagonal blocks
 * during the NSGS iterations.
 *
 * For other matrix types, this function does nothing.
 *
 * \param[in,out] M The matrix to extract diagonal blocks from
 * \param[in] dim Dimension per block
 * \return Pointer to the original M->matrix1 (for later restoration), or NULL if no extraction
 * was done
 */
static inline void* nsgs_generic_extract_diagonal_blocks(NumericsMatrix* M, int dim) {
  if (!M) return NULL;

  if (M->storageType == NM_SPARSE) {
    void* original_matrix1 = M->matrix1;
    if (M->matrix1) {
      numerics_printf(
          "nsgs_generic_extract_diagonal_blocks: Warning - M->matrix1 already exists, "
          "overwriting");
    }
    M->matrix1 = NM_extract_diagonal_blocks(M, dim);
    return original_matrix1;
  }

  return NULL;
}

/** Free extracted diagonal blocks and restore original matrix1 pointer
 *
 * Cleans up the memory allocated by nsgs_generic_extract_diagonal_blocks
 * and restores the original M->matrix1 pointer.
 *
 * \param[in,out] M The matrix to clean up
 * \param[in] original_matrix1 The original M->matrix1 pointer (from
 * nsgs_generic_extract_diagonal_blocks)
 */
static inline void nsgs_generic_free_diagonal_blocks(NumericsMatrix* M,
                                                     void* original_matrix1) {
  if (!M) return;

  if (M->storageType == NM_SPARSE && M->matrix1) {
    /* Free the extracted diagonal blocks */
    SBM_clear_block(M->matrix1);
    SBM_clear(M->matrix1);
    free(M->matrix1);
    /* Restore original matrix1 */
    M->matrix1 = (SparseBlockStructuredMatrix*)original_matrix1;
  }
}

/** Determine convergence based on error and tolerance
 *
 * \param error Current error
 * \param tolerance Convergence tolerance
 * \param iter Current iteration
 * \param options Solver options
 * \return 1 if not converged, 0 if converged
 */
static inline int nsgs_determine_convergence(double error, double tolerance, int iter,
                                             SolverOptions* options) {
  (void)iter; /* May be used for adaptive tolerance */
  (void)options;
  return (error > tolerance) ? 1 : 0;
}

/** Check if too many blocks are frozen and reset if necessary
 *
 * If nb_blocks - 1 or more blocks are frozen, reset all to unfreeze
 * at least one block for the next iteration.
 *
 * \param freeze_blocks Array of freeze counters
 * \param nb_blocks Number of blocks
 * \param use_freezing Whether freezing is enabled
 * \return Number of frozen blocks
 */
static inline unsigned int nsgs_check_freezing_reset(unsigned int* freeze_blocks,
                                                     unsigned int nb_blocks,
                                                     int use_freezing) {
  if (!use_freezing || !freeze_blocks) return 0;

  unsigned int nb_frozen = 0;
  for (unsigned int i = 0; i < nb_blocks; ++i) {
    if (freeze_blocks[i] > 0) nb_frozen++;
  }

  /* If too many blocks frozen, reset all */
  if (nb_frozen >= nb_blocks - 1) {
    for (unsigned int i = 0; i < nb_blocks; ++i) {
      freeze_blocks[i] = 0;
    }
    nb_frozen = 0;
  }

  return nb_frozen;
}

/** Default criteria for freezing a block (matches fc3d_nsgs)
 *
 * A block is frozen if either:
 * 1. Relative convergence: incremental_error <= tolerance^2 * 100^2 * ||z_local||^2
 * 2. Small solution: ||z_local||^2 <= ||z_global||^2 / (nb_blocks^2 * 1000)
 *
 * And we are past iteration 10.
 *
 * \param incremental_error_2 Squared incremental error for this block (||z_new - z_old||^2)
 * \param squared_norm_local Squared norm of local solution (||z_local||^2)
 * \param tolerance Global tolerance
 * \param norm_z_global Norm of global solution
 * \param nb_blocks Number of blocks
 * \param iter Current iteration
 * \return 1 if block should be frozen, 0 otherwise
 */
static inline int nsgs_default_should_freeze(double incremental_error_2,
                                             double squared_norm_local, double tolerance,
                                             double norm_z_global, unsigned int nb_blocks,
                                             int iter) {
  if (iter < 10) return 0;
  if (nb_blocks == 0) return 0;

  /* Criteria matching fc3d_nsgs exactly */
  double tmp_criteria1 = tolerance * tolerance * 100.0 * 100.0;
  double tmp_criteria2 = norm_z_global * norm_z_global / (nb_blocks * nb_blocks * 1000.0);

  int relative_convergence_criteria =
      incremental_error_2 <= tmp_criteria1 * squared_norm_local;
  int small_reaction_criteria = squared_norm_local <= tmp_criteria2;

  return (relative_convergence_criteria || small_reaction_criteria) ? 1 : 0;
}

/** Update freeze counters for a block
 *
 * \param freeze_blocks Array of freeze counters
 * \param block Block index
 * \param freezing_iter Number of iterations to freeze
 * \param should_freeze Whether this block should be frozen
 */
static inline void nsgs_update_freeze_counter(unsigned int* freeze_blocks, unsigned int block,
                                              int freezing_iter, int should_freeze) {
  if (!freeze_blocks) return;

  if (should_freeze) {
    freeze_blocks[block] = freezing_iter;
  }
}

/** Shuffle blocks array
 *
 * \param sblocks Array of block indices to shuffle
 * \param nb_blocks Number of blocks
 * \param shuffle_mode Shuffle mode (0=none, 1=once, 2=each loop)
 * \param iter Current iteration
 */
static inline void nsgs_shuffle_blocks(unsigned int* sblocks, unsigned int nb_blocks,
                                       int shuffle_mode, int iter) {
  if (!sblocks) return;
  if (shuffle_mode == 0) return;

  /* Shuffle each loop, or first iteration only */
  if (shuffle_mode == 2 || (shuffle_mode == 1 && iter == 1)) {
    /* Note: uint_shuffle is from NumericsArrays.h */
    extern void uint_shuffle(unsigned int*, unsigned int);
    uint_shuffle(sblocks, nb_blocks);
  }
}

/** Count frozen blocks and return percentage
 *
 * \param[in] freeze_blocks Array of freeze counters
 * \param[in] nb_blocks Number of blocks
 * \param[in] use_freezing Whether freezing is enabled
 * \return Percentage of frozen blocks (0-100)
 */
static inline double nsgs_count_frozen_percent(unsigned int* freeze_blocks,
                                               unsigned int nb_blocks, int use_freezing) {
  if (!use_freezing || !freeze_blocks || nb_blocks == 0) return 0.0;

  unsigned int nb_frozen = 0;
  for (unsigned int i = 0; i < nb_blocks; ++i) {
    if (freeze_blocks[i] > 0) nb_frozen++;
  }
  return ((double)nb_frozen / nb_blocks) * 100.0;
}

/** Compute full error and check convergence with user tolerance
 *
 * This function computes the full error and determines if the solution
 * has converged based on the user-specified tolerance.
 *
 * \param[in] problem Problem data
 * \param[in] var_z Global solution vector
 * \param[in] var_x Global dual vector
 * \param[in] options Solver options
 * \param[in] toolkit Local problem toolkit
 * \param[in] user_tolerance User-specified tolerance
 * \param[out] full_error Computed full error
 * \return 0 if converged (full_error <= user_tolerance), 1 otherwise
 */
static inline int nsgs_check_full_error_convergence(void* problem, double* var_z,
                                                    double* var_x, SolverOptions* options,
                                                    NSGSLocalToolkit* toolkit,
                                                    double user_tolerance,
                                                    double* full_error) {
  if (!toolkit->compute_error) {
    *full_error = 0.0;
    return 0;
  }

  *full_error = toolkit->compute_error(problem, var_z, var_x, options);

  return (*full_error <= user_tolerance) ? 0 : 1;
}

/** Early tolerance adaptation at iteration 10
 *
 * Compute full error and adapt working tolerance based on incr/full ratio.
 *
 * \param[in] problem Problem data
 * \param[in] var_z Global solution vector
 * \param[in] var_x Global dual vector
 * \param[in] options Solver options
 * \param[in] toolkit Local problem toolkit
 * \param[in] incremental_error Current incremental error
 * \param[in,out] tolerance Working tolerance (may be adapted)
 * \param[out] full_error Computed full error
 * \param[in] iter Current iteration
 */
static inline void nsgs_early_tolerance_adapt(void* problem, double* var_z, double* var_x,
                                              SolverOptions* options,
                                              NSGSLocalToolkit* toolkit,
                                              double incremental_error, double* tolerance,
                                              double* full_error, int iter) {
  if (iter != 10 || !toolkit->compute_error) return;

  nsgs_check_full_error_convergence(problem, var_z, var_x, options, toolkit,
                                    toolkit->user_tolerance, full_error);

  /* Adapt tolerance based on the ratio at iteration 10 (use ratio directly) */
  if (incremental_error > 0 && *full_error > 0) {
    double error_ratio = incremental_error / *full_error;
    double new_tolerance = error_ratio * toolkit->user_tolerance;

    if (toolkit->verbose > 0) {
      printf(
          "Iter %d: Early adapt - Incr: %.4e, Full: %.4e, Ratio: %.4f, "
          "Current tol: %.4e, New tol: %.4e %s\n",
          iter, incremental_error, *full_error, error_ratio, *tolerance, new_tolerance,
          (new_tolerance != *tolerance) ? "[APPLYING]" : "[NO CHANGE]");
    }

    /* Always apply the new tolerance (as long as it's positive) */
    if (new_tolerance > 0) {
      *tolerance = new_tolerance;
    }
  }
}

/** Print iteration statistics (verbose output)
 *
 * \param[in] iter Current iteration
 * \param[in] incremental_error incremental error value
 * \param[in] full_error Full error value
 * \param[in] tolerance Working tolerance
 * \param[in] user_tolerance User tolerance
 * \param[in] frozen_percent Percentage of frozen blocks
 * \param[in] hasNotConverged Convergence flag
 * \param[in] verbose Verbosity level
 */
static inline void nsgs_print_iteration_stats(int iter, double incremental_error,
                                              double full_error, double tolerance,
                                              double user_tolerance, double frozen_percent,
                                              int hasNotConverged, int verbose) {
  if (verbose <= 0 || iter % verbose != 0) return;

  /* Print header on first iteration */
  if (iter == 1) {
    printf("\n%-8s %-12s %-12s %-12s %-12s %-8s %-8s\n", "Iter", "IncrError", "FullError",
           "WorkTol", "UserTol", "Frozen%", "Status");
    printf("%-8s %-12s %-12s %-12s %-12s %-8s %-8s\n", "--------", "------------",
           "------------", "------------", "------------", "--------", "--------");
  }

  const char* status = hasNotConverged ? "..." : "CONV";

  printf("%-8d %-12.4e %-12.4e %-12.4e %-12.4e %-7.1f%% %-8s\n", iter, incremental_error,
         full_error, tolerance, user_tolerance, frozen_percent, status);
}

/** Adapt working tolerance based on incremental and full error
 *
 * When incremental error indicates convergence but full error does not,
 * adapt the working tolerance to try to achieve full convergence.
 *
 * Strategy:
 * - If incremental_error is very small (< 1e-15): tighten local solver
 * - Otherwise: adapt working tolerance based on error ratio
 *
 * \param[in] incremental_error Current incremental error (squared norm of change)
 * \param[in] full_error Current full error
 * \param[in] user_tolerance User-specified tolerance
 * \param[in,out] working_tolerance Current working tolerance (may be adapted)
 * \param[in] localsolver_options Internal solver options (may be tightened)
 * \param[in] iter Current iteration
 * \param[in] verbose Verbosity level
 * \return 1 if we should continue iterating, 0 if converged
 */
static inline int nsgs_adapt_tolerance(double incremental_error, double full_error,
                                       double user_tolerance, double* working_tolerance,
                                       SolverOptions* localsolver_options, int iter,
                                       int verbose) {
  /* If full error is within user tolerance, we're converged */
  if (full_error <= user_tolerance) {
    if (verbose > 0) {
      printf("Iter %d: CONVERGED - Full error (%.2e) <= user tolerance (%.2e)\n", iter,
             full_error, user_tolerance);
    }
    return 0; /* Converged */
  }

  /* Full error not converged - need to adapt */
  if (verbose > 0) {
    printf(
        "Iter %d: Incr error (%.2e) < work_tol (%.2e), but full error (%.2e) > user_tol "
        "(%.2e)\n",
        iter, incremental_error, *working_tolerance, full_error, user_tolerance);
  }

  /* If incremental error is very small, the NSGS loop is not improving the solution.
   * We need to tighten the local solver tolerance instead of adapting working tolerance. */
  if (incremental_error < 1e-15 || incremental_error < DBL_EPSILON * 100.0) {
    if (localsolver_options) {
      double current_local_tol = localsolver_options->dparam[SICONOS_DPARAM_TOL];
      double new_local_tol = fmax(current_local_tol / 100.0, DBL_EPSILON * 1e-6);
      localsolver_options->dparam[SICONOS_DPARAM_TOL] = new_local_tol;
      if (verbose > 0) {
        printf(
            "Iter %d: Incr error very small (%.2e), tightening local solver: %.2e -> %.2e\n",
            iter, incremental_error, current_local_tol, new_local_tol);
      }
    } else if (verbose > 0) {
      printf("Iter %d: Incr error very small (%.2e) but no local solver options to tighten\n",
             iter, incremental_error);
    }
    /* Don't change working tolerance when incr error is negligible */
    return 1; /* Not converged, continue iterating */
  }

  /* Normal case: adapt working tolerance based on error ratio */
  if (full_error > 0 && incremental_error > 0) {
    /* The ratio tells us how much tighter we need to go */
    double error_ratio = incremental_error / full_error;
    double new_tolerance = error_ratio * user_tolerance;

    /* Sanity checks: don't go below machine epsilon, don't increase tolerance */
    double min_tolerance = user_tolerance * 1e-6;
    if (new_tolerance < min_tolerance) {
      if (verbose > 1) {
        printf("Iter %d: Capping tolerance at %.2e (computed: %.2e)\n", iter, min_tolerance,
               new_tolerance);
      }
      new_tolerance = min_tolerance;
    }

    if (new_tolerance > 0 && new_tolerance < *working_tolerance) {
      if (verbose > 0) {
        printf("Iter %d: Adapting working tolerance: %.4e -> %.4e (ratio=%.4f)\n", iter,
               *working_tolerance, new_tolerance, error_ratio);
      }
      *working_tolerance = new_tolerance;
    } else if (verbose > 1) {
      printf("Iter %d: Skipping tolerance adaptation (new=%.2e not < current=%.2e)\n", iter,
             new_tolerance, *working_tolerance);
    }
  }

  return 1; /* Not converged, continue iterating */
}
/** Check convergence when incremental error has converged
 *
 * If incremental error <= tolerance, compute full error and check against user tolerance.
 * Adapt tolerance if incremental converged but full did not.
 *
 * \param[in] problem Problem data
 * \param[in] var_z Global solution vector
 * \param[in] var_x Global dual vector
 * \param[in] options Solver options
 * \param[in] toolkit Local problem toolkit
 * \param[in] incremental_error Current incremental error
 * \param[in,out] tolerance Working tolerance (may be adapted)
 * \param[in,out] full_error Computed full error (output)
 * \param[in] computed_full_error Whether full error was already computed
 * \param[in] localsolver_options Internal solver options
 * \param[in] iter Current iteration
 * \return 0 if converged, 1 if not converged
 */
static inline int nsgs_check_convergence_with_full_error(
    void* problem, double* var_z, double* var_x, SolverOptions* options,
    NSGSLocalToolkit* toolkit, double incremental_error, double* tolerance, double* full_error,
    int computed_full_error, SolverOptions* localsolver_options, int iter) {
  /* Compute full error if not already done */
  if (!computed_full_error) {
    nsgs_check_full_error_convergence(problem, var_z, var_x, options, toolkit,
                                      toolkit->user_tolerance, full_error);
  }

  /* Check if full error is within user tolerance */
  if (*full_error <= toolkit->user_tolerance) {
    /* Both incremental and full errors converged */
    if (toolkit->verbose > 0) {
      printf(
          "Iter %d: CONVERGED - Incr (%.2e) <= tol (%.2e), Full (%.2e) <= user_tol (%.2e)\n",
          iter, incremental_error, *tolerance, *full_error, toolkit->user_tolerance);
    }
    return 0; /* Converged */
  }

  /* Incremental converged but full did not - adapt tolerance and continue */
  return nsgs_adapt_tolerance(incremental_error, *full_error, toolkit->user_tolerance,
                              tolerance, localsolver_options, iter, toolkit->verbose);
}
/** Generic NSGS solver - SINGLE LOOP VERSION (Optimized)
 *
 * This version uses a single loop with inline freezing and error computation,
 * matching the performance characteristics of the original fc3d_nsgs.
 *
 * Key optimizations:
 * - Single pass over blocks (no separate freezing stage)
 * - Per-block error accumulation (no diff_vec needed)
 * - Inline small vector operations (no BLAS overhead for dim <= 5)
 * - Immediate freezing decision (affects current iteration)
 *
 * \param[in] problem Problem data (problem-specific structure)
 * \param[in,out] var_z  solution vector (primal, e.g., reaction)
 * \param[in,out] var_x  dual vector (e.g., velocity, can be NULL)
 * \param[out] info Solver status (0=success, 1=failure)
 * \param[in] options Solver options
 * \param[in] toolkit Local problem toolkit with callbacks
 * \param[in] problem_data Common problem data
 */
static inline void nsgs_solve(void* problem, double* var_z, double* var_x, int* info,
                              SolverOptions* options, NSGSLocalToolkit* toolkit,
                              NSGSProblemData* problem_data) {
  /* Get solver parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  unsigned int nb_blocks = problem_data->nb_blocks;
  int dim = toolkit->dimension;

  // double norm_q = cblas_dnrm2(nb_blocks * dim, problem_data->q, 1);
  double norm_z = 1e24;
  double prev_norm_z = 1e24; /* For freezing criteria (like original) */

  /* Iteration variables */
  int iter = 0;
  double error = 1.0;
  int hasNotConverged = 1;

  /* Allocate working buffers - minimal for single-loop version */
  double* var_z_local = (double*)malloc(dim * sizeof(double));
  double* block_errors = NULL; /* Per-block squared errors for freezing */
  if (toolkit->use_incremental_error || toolkit->use_freezing) {
    block_errors = (double*)malloc(nb_blocks * sizeof(double));
  }

  /* Get local problem: use alloc_local callback if provided, otherwise use
   * toolkit.localproblem */
  void* local_problem = NULL;
  /* Allocate block ordering arrays */
  unsigned int* sblocks = NULL;
  unsigned int* freeze_blocks = NULL;
  SolverOptions* localsolver_options = NULL;

  if (toolkit->alloc_local) {
    local_problem = toolkit->alloc_local(problem);
    if (!local_problem) {
      *info = numerics_error("nsgs_solve", "Failed to allocate local problem");
      goto nsgs_cleanup;
    }
  } else if (toolkit->localproblem) {
    /* Use pre-allocated local problem from toolkit */
    local_problem = toolkit->localproblem;
  }

  /* Allocate shuffled blocks array if shuffling is enabled */
  if (toolkit->use_shuffling) {
    if (toolkit->alloc_shuffled) {
      sblocks = toolkit->alloc_shuffled(problem, options);
    } else {
      /* Default allocation: allocate and initialize with 0, 1, 2, ... */
      sblocks = (unsigned int*)malloc(nb_blocks * sizeof(unsigned int));
      for (unsigned int i = 0; i < nb_blocks; ++i) {
        sblocks[i] = i;
      }
    }
  }

  /* Allocate freeze blocks array if freezing is enabled */
  if (toolkit->use_freezing) {
    if (toolkit->alloc_freezing) {
      freeze_blocks = toolkit->alloc_freezing(problem, options);
    } else {
      /* Default allocation: allocate and initialize to zero */
      freeze_blocks = (unsigned int*)calloc(nb_blocks, sizeof(unsigned int));
    }
  }

  /* Get internal solver options */
  double local_tolerance_save = 0.0;
  if (options->numberOfInternalSolvers > 0) {
    localsolver_options = options->internalSolvers[0];
    local_tolerance_save = localsolver_options->dparam[SICONOS_DPARAM_TOL];
  }

  /* Check for early exit - if *info == 0 on entry, problem is already solved (trivial case) */
  if (*info == 0) {
    goto nsgs_cleanup;
  }

  /*****  Main NSGS Iterations - SINGLE LOOP VERSION *****/
  while ((iter < itermax) && hasNotConverged) {
    ++iter;
    double incremental_error_sum = 0.0;

    /* Set tolerance for internal solver */
    if (toolkit->set_local_tol && localsolver_options) {
      toolkit->set_local_tol(problem, options, localsolver_options, error);
    }

    /* Check freezing state - reset if too many blocks frozen */
    nsgs_check_freezing_reset(freeze_blocks, nb_blocks, toolkit->use_freezing);

    /* Shuffle blocks if needed */
    nsgs_shuffle_blocks(sblocks, nb_blocks, toolkit->use_shuffling ? 2 : 0, iter);

    /* Pre-compute freezing criteria constants (like original fc3d_nsgs) */
    double tmp_criteria1 = tolerance * tolerance  / (nb_blocks * nb_blocks * 1000.0);
    double tmp_criteria2 = (prev_norm_z > 0.0)
                               ? (prev_norm_z * prev_norm_z / (nb_blocks * nb_blocks * 1000.0))
                               : 0.0;

    /* =====================================================================
     * SINGLE LOOP: Process all blocks with inline error and freezing
     * ===================================================================== */
    for (unsigned int i = 0; i < nb_blocks; ++i) {
      unsigned int block;

      /* Determine block index (shuffled or sequential) */
      if (sblocks && toolkit->use_shuffling) {
        block = sblocks[i];
      } else {
        block = i;
      }

      /* Skip frozen blocks */
      if (freeze_blocks && toolkit->use_freezing && freeze_blocks[block] > 0) {
        freeze_blocks[block]--;
        continue;
      }

      /* Copy current solution to local buffer - inline for small dim */
      if (dim <= 5) {
        for (int d = 0; d < dim; ++d) {
          var_z_local[d] = var_z[block * dim + d];
        }
      } else {
        cblas_dcopy(dim, &var_z[block * dim], 1, var_z_local, 1);
      }

      /* Update local problem */
      if (toolkit->update_local_problem) {
        toolkit->update_local_problem(block, problem, local_problem, var_z, options);
      }

      /* Solve local problem */
      if (toolkit->solve_local) {
        toolkit->solve_local(local_problem, localsolver_options, var_z_local, NULL);
      }

      /* Perform relaxation if enabled */
      if (toolkit->use_relaxation && toolkit->relaxation) {
        toolkit->relaxation(var_z_local, &var_z[block * dim], toolkit->omega);
      }

      /* Compute incremental error contribution and squared norm - inline */
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
          block_errors[block] = incr_err_2; /* Store for potential reuse */
        }
      }

      /* Check if block should be frozen - IMMEDIATE (like original fc3d_nsgs) */
      if (freeze_blocks && toolkit->use_freezing && toolkit->freezing_iter > 0 && iter >= 10) {
        int relative_conv = incr_err_2 <= tmp_criteria1 * sq_norm_local;
        int small_sol = sq_norm_local <= tmp_criteria2;
        if (relative_conv || small_sol) {
          freeze_blocks[block] = toolkit->freezing_iter;
        }
      }

      /* Accept local solution into global vector - inline for small dim */
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
    }

    /* =====================================================================
     * Compute global error and check convergence
     * ===================================================================== */
    prev_norm_z = norm_z; /* Save for next iteration's freezing criteria */
    norm_z = cblas_dnrm2(nb_blocks * dim, var_z, 1);

    double incremental_error = 0.0;
    if (toolkit->use_incremental_error) {
      incremental_error = sqrt(incremental_error_sum);
      if (norm_z > 0.0) {
        incremental_error /= norm_z;
      }
    }

    /* Initialize error tracking */
    hasNotConverged = 1;
    error = incremental_error;
    double full_error = 0.0;

    /* Note: Original fc3d_nsgs does NOT do early tolerance adaptation at iter 10.
     * Tolerance is only adapted when incremental error converges but full error doesn't.
     * This matches the behavior in determine_convergence_with_full_final(). */

    /* Check convergence based on incremental error */
    if (incremental_error <= tolerance) {
      /* Incremental converged - check full error and adapt tolerance if needed
       * This matches fc3d_nsgs behavior: only adapt when incr converges but full doesn't */
      int full_already_computed = 0;
      hasNotConverged = nsgs_check_convergence_with_full_error(
          problem, var_z, var_x, options, toolkit, incremental_error, &tolerance, &full_error,
          full_already_computed, localsolver_options, iter);
      error = full_error;
    } else {
      /* Incremental not converged - standard check */
      if (toolkit->check_convergence) {
        hasNotConverged =
            toolkit->check_convergence(incremental_error, tolerance, iter, options);
      } else {
        hasNotConverged =
            nsgs_determine_convergence(incremental_error, tolerance, iter, options);
      }
    }

    /* Print iteration statistics */
    double frozen_percent =
        nsgs_count_frozen_percent(freeze_blocks, nb_blocks, toolkit->use_freezing);

    nsgs_print_iteration_stats(iter, incremental_error, full_error, tolerance,
                               toolkit->user_tolerance, frozen_percent, hasNotConverged,
                               toolkit->verbose);

    /* Statistics callback */
    if (toolkit->stats_callback) {
      toolkit->stats_callback(problem, options, var_z, var_x, error);
    }
  }

  /* Print newline after iterations complete if verbose was enabled */
  if (toolkit->verbose > 0 && iter > 0) {
    printf("\n");
  }

  /* Final error verification: always check full error if compute_error is available */
  if (toolkit->compute_error) {
    double final_full_error;
    int final_converged = !nsgs_check_full_error_convergence(
        problem, var_z, var_x, options, toolkit, toolkit->user_tolerance, &final_full_error);

    if (toolkit->verbose > 0) {
      printf("Final check: Full error = %.6e, User tolerance = %.6e, Converged = %s\n",
             final_full_error, toolkit->user_tolerance, final_converged ? "YES" : "NO");
    }

    /* Final result based on full error vs user tolerance */
    error = final_full_error;
    *info = final_converged ? 0 : 1;
    hasNotConverged = final_converged ? 0 : 1;
  }

  /* Set output info */
  if (error <= tolerance) {
    *info = 0;
  } else {
    *info = 1;
  }

  /* Store iteration info */
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  dparam[SICONOS_DPARAM_RESIDU] = error;

nsgs_cleanup:
  /* Restore local solver tolerance */
  if (localsolver_options) {
    localsolver_options->dparam[SICONOS_DPARAM_TOL] = local_tolerance_save;
  }

  /* Cleanup */
  free(var_z_local);
  if (block_errors) free(block_errors);
  if (sblocks) free(sblocks);
  if (freeze_blocks) free(freeze_blocks);
  /* Note: local_problem cleanup is caller's responsibility if alloc_local was used */
}

/** Generic NSGS solver - 4 STAGE VERSION (Original Design)
 *
 * This version separates the computation into 4 stages:
 * 1. Local solutions for all blocks
 * 2. Error computation
 * 3. Freezing status update
 * 4. Statistics printing
 *
 * This design is more modular but has higher overhead.
 *
 * \param[in] problem Problem data (problem-specific structure)
 * \param[in,out] var_z  solution vector (primal, e.g., reaction)
 * \param[in,out] var_x  dual vector (e.g., velocity, can be NULL)
 * \param[out] info Solver status (0=success, 1=failure)
 * \param[in] options Solver options
 * \param[in] toolkit Local problem toolkit with callbacks
 * \param[in] problem_data Common problem data
 */
static inline void nsgs_solve_4stage(void* problem, double* var_z, double* var_x, int* info,
                                     SolverOptions* options, NSGSLocalToolkit* toolkit,
                                     NSGSProblemData* problem_data) {
  /* Get solver parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  unsigned int nb_blocks = problem_data->nb_blocks;
  int dim = toolkit->dimension;

  // double norm_q = cblas_dnrm2(nb_blocks * dim, problem_data->q, 1);
  double norm_z = 1e24;

  /* Iteration variables */
  int iter = 0;
  double error = 1.0;
  int hasNotConverged = 1;

  /* Allocate working buffers */
  double* var_z_local = (double*)malloc(dim * sizeof(double));
  double* var_z_old = NULL;
  double* diff_vec = NULL; /* For BLAS-based incremental error computation */
  if (toolkit->use_incremental_error || toolkit->use_freezing) {
    var_z_old = (double*)malloc(nb_blocks * dim * sizeof(double));
    if (toolkit->use_incremental_error) {
      diff_vec = (double*)malloc(nb_blocks * dim * sizeof(double));
    }
  }

  /* Get local problem: use alloc_local callback if provided, otherwise use
   * toolkit.localproblem */
  void* local_problem = NULL;
  if (toolkit->alloc_local) {
    local_problem = toolkit->alloc_local(problem);
    if (!local_problem) {
      *info = numerics_error("nsgs_solve", "Failed to allocate local problem");
      goto nsgs_cleanup;
    }
  } else if (toolkit->localproblem) {
    /* Use pre-allocated local problem from toolkit */
    local_problem = toolkit->localproblem;
  }

  /* Allocate block ordering arrays */
  unsigned int* sblocks = NULL;
  unsigned int* freeze_blocks = NULL;

  /* Allocate shuffled blocks array if shuffling is enabled */
  if (toolkit->use_shuffling) {
    if (toolkit->alloc_shuffled) {
      sblocks = toolkit->alloc_shuffled(problem, options);
    } else {
      /* Default allocation: allocate and initialize with 0, 1, 2, ... */
      sblocks = (unsigned int*)malloc(nb_blocks * sizeof(unsigned int));
      for (unsigned int i = 0; i < nb_blocks; ++i) {
        sblocks[i] = i;
      }
    }
  }

  /* Allocate freeze blocks array if freezing is enabled */
  if (toolkit->use_freezing) {
    if (toolkit->alloc_freezing) {
      freeze_blocks = toolkit->alloc_freezing(problem, options);
    } else {
      /* Default allocation: allocate and initialize to zero */
      freeze_blocks = (unsigned int*)calloc(nb_blocks, sizeof(unsigned int));
    }
  }

  /* Get internal solver options */
  SolverOptions* localsolver_options = NULL;
  if (options->numberOfInternalSolvers > 0) {
    localsolver_options = options->internalSolvers[0];
  }

  /* Check for early exit - if *info == 0 on entry, problem is already solved (trivial case) */
  if (*info == 0) {
    goto nsgs_cleanup;
  }

  /*****  Main NSGS Iterations *****/
  while ((iter < itermax) && hasNotConverged) {
    ++iter;

    /* Set tolerance for internal solver */
    if (toolkit->set_local_tol && localsolver_options) {
      toolkit->set_local_tol(problem, options, localsolver_options, error);
    }

    /* Check freezing state - reset if too many blocks frozen */
    nsgs_check_freezing_reset(freeze_blocks, nb_blocks, toolkit->use_freezing);

    /* Shuffle blocks if needed */
    nsgs_shuffle_blocks(sblocks, nb_blocks, toolkit->use_shuffling ? 2 : 0, iter);

    /* =====================================================================
     * STAGE 1: Compute local solutions for all blocks
     * ===================================================================== */
    for (unsigned int i = 0; i < nb_blocks; ++i) {
      unsigned int block;

      /* Determine block index (shuffled or sequential) */
      if (sblocks && toolkit->use_shuffling) {
        block = sblocks[i];
      } else {
        block = i;
      }

      /* Skip frozen blocks */
      if (freeze_blocks && toolkit->use_freezing && freeze_blocks[block] > 0) {
        freeze_blocks[block]--;
        continue;
      }

      /* Store old solution for error computation and freezing using BLAS */
      if (var_z_old) {
        cblas_dcopy(dim, &var_z[block * dim], 1, &var_z_old[block * dim], 1);
      }

      /* Copy current solution to local buffer using BLAS */
      cblas_dcopy(dim, &var_z[block * dim], 1, var_z_local, 1);

      /* Update local problem */
      if (toolkit->update_local_problem) {
        toolkit->update_local_problem(block, problem, local_problem, var_z, options);
      }

      /* Solve local problem */
      if (toolkit->solve_local) {
        toolkit->solve_local(local_problem, localsolver_options, var_z_local, NULL);
      }

      /* Perform relaxation if enabled */
      if (toolkit->use_relaxation && toolkit->relaxation) {
        toolkit->relaxation(var_z_local, &var_z[block * dim], toolkit->omega);
      }

      /* Accept local solution into global vector */
      if (toolkit->accept_local) {
        toolkit->accept_local(local_problem, localsolver_options, block, iter, var_z,
                              var_z_local);
      } else {
        /* Default: copy local solution to global using BLAS */
        cblas_dcopy(dim, var_z_local, 1, &var_z[block * dim], 1);
      }
    }

    /* =====================================================================
     * STAGE 2: Compute error(s)
     * ===================================================================== */
    double incremental_error = 0.0;
    double incr_err_sq = 0.0; /* Squared incremental error for freezing */
    norm_z = cblas_dnrm2(nb_blocks * dim, var_z, 1);

    if (toolkit->use_incremental_error && var_z_old) {
      /* Optimized incremental error computation using BLAS:
       * 1. diff_vec = var_z - var_z_old  (using daxpy: y = -1*old + z)
       * 2. incr_err_sq = diff_vec^T * diff_vec  (using ddot)
       * 3. incremental_error = sqrt(incr_err_sq) / norm_z
       */
      if (diff_vec) {
        /* diff_vec = var_z */
        cblas_dcopy(nb_blocks * dim, var_z, 1, diff_vec, 1);
        /* diff_vec = diff_vec - var_z_old = var_z - var_z_old */
        cblas_daxpy(nb_blocks * dim, -1.0, var_z_old, 1, diff_vec, 1);
        /* incr_err_sq = diff_vec^T * diff_vec */
        incr_err_sq = cblas_ddot(nb_blocks * dim, diff_vec, 1, diff_vec, 1);
        incremental_error = sqrt(incr_err_sq);
        if (norm_z > 0.0) {
          incremental_error /= norm_z;
        }
      }
    }

    /* Initialize error tracking */
    hasNotConverged = 1;
    error = incremental_error;
    double full_error = 0.0;
    int computed_full_error = 0;

    /* At iteration 10, compute full error to adapt tolerance early */
    nsgs_early_tolerance_adapt(problem, var_z, var_x, options, toolkit, incremental_error,
                               &tolerance, &full_error, iter);
    if (iter == 10) computed_full_error = 1;

    /* Check convergence based on incremental error */
    if (incremental_error <= tolerance) {
      /* Incremental converged - check full error and adapt if needed */
      hasNotConverged = nsgs_check_convergence_with_full_error(
          problem, var_z, var_x, options, toolkit, incremental_error, &tolerance, &full_error,
          computed_full_error, localsolver_options, iter);
      error = full_error;
    } else {
      /* Incremental not converged - standard check */
      if (toolkit->check_convergence) {
        hasNotConverged =
            toolkit->check_convergence(incremental_error, tolerance, iter, options);
      } else {
        hasNotConverged =
            nsgs_determine_convergence(incremental_error, tolerance, iter, options);
      }
    }

    /* =====================================================================
     * STAGE 3: Update freezing status
     * ===================================================================== */
    if (freeze_blocks && toolkit->use_freezing && toolkit->freezing_iter > 0 && var_z_old) {
      for (unsigned int block = 0; block < nb_blocks; ++block) {
        /* Compute per-block incremental error using diff_vec if available */
        double incr_err_2 = 0.0;
        if (diff_vec) {
          /* Use pre-computed diff_vec from Stage 2 */
          incr_err_2 = cblas_ddot(dim, &diff_vec[block * dim], 1, &diff_vec[block * dim], 1);
        } else {
          /* Fallback: compute diff on the fly */
          for (int d = 0; d < dim; ++d) {
            double diff = var_z[block * dim + d] - var_z_old[block * dim + d];
            incr_err_2 += diff * diff;
          }
        }

        /* Compute squared norm of local solution using BLAS or callback */
        double sq_norm_local = 0.0;
        if (toolkit->squared_norm) {
          sq_norm_local = toolkit->squared_norm(&var_z[block * dim]);
        } else {
          sq_norm_local = cblas_ddot(dim, &var_z[block * dim], 1, &var_z[block * dim], 1);
        }

        int should_freeze = 0;
        if (toolkit->should_freeze) {
          should_freeze = toolkit->should_freeze(incr_err_2, sq_norm_local, tolerance, norm_z,
                                                 nb_blocks, iter);
        } else {
          should_freeze = nsgs_default_should_freeze(incr_err_2, sq_norm_local, tolerance,
                                                     norm_z, nb_blocks, iter);
        }

        if (should_freeze) {
          freeze_blocks[block] = toolkit->freezing_iter;
        }
      }
    }

    /* =====================================================================
     * STAGE 4: Print iteration statistics
     * ===================================================================== */
    double frozen_percent =
        nsgs_count_frozen_percent(freeze_blocks, nb_blocks, toolkit->use_freezing);

    nsgs_print_iteration_stats(iter, incremental_error, full_error, tolerance,
                               toolkit->user_tolerance, frozen_percent, hasNotConverged,
                               toolkit->verbose);

    /* Statistics callback */
    if (toolkit->stats_callback) {
      toolkit->stats_callback(problem, options, var_z, var_x, error);
    }
  }

  /* Print newline after iterations complete if verbose was enabled */
  if (toolkit->verbose > 0 && iter > 0) {
    printf("\n");
  }

  /* Final error verification: always check full error if compute_error is available */
  if (toolkit->compute_error) {
    double final_full_error;
    int final_converged = !nsgs_check_full_error_convergence(
        problem, var_z, var_x, options, toolkit, toolkit->user_tolerance, &final_full_error);

    if (toolkit->verbose > 0) {
      printf("Final check: Full error = %.6e, User tolerance = %.6e, Converged = %s\n",
             final_full_error, toolkit->user_tolerance, final_converged ? "YES" : "NO");
    }

    /* Final result based on full error vs user tolerance */
    error = final_full_error;
    *info = final_converged ? 0 : 1;
    hasNotConverged = final_converged ? 0 : 1;
  }

  /* Set output info */
  if (error <= tolerance) {
    *info = 0;
  } else {
    *info = 1;
  }

  /* Store iteration info */
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  dparam[SICONOS_DPARAM_RESIDU] = error;

nsgs_cleanup:
  /* Cleanup */
  free(var_z_local);
  if (var_z_old) free(var_z_old);
  if (diff_vec) free(diff_vec);
  if (sblocks) free(sblocks);
  if (freeze_blocks) free(freeze_blocks);
  /* Note: local_problem cleanup is caller's responsibility if alloc_local was used */
}

#endif /* NSGS_GENERIC_H */
