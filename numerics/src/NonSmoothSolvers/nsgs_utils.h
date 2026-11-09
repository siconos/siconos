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

#ifndef NSGS_UTILS_H
#define NSGS_UTILS_H

/*!\file nsgs_utils.h
 * \brief Problem-independent helpers shared by NSGS implementations.
 *
 * These small inline utilities were factored out of nsgs_generic.h so that
 * problem-specific NSGS solvers can reuse them without importing the full
 * generic solver framework.
 */

#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs
#include <stdio.h>   // for printf
#include <stdlib.h>  // for size_t

#include "NumericsArrays.h"  // for uint_shuffle
#include "SolverOptions.h"   // for SolverOptions, SICONOS_DPARAM_TOL

/** Relaxation: z_new = omega * z_new + (1 - omega) * z_old.
 *
 * \param[in,out] z_new new local solution, overwritten with relaxed value
 * \param[in] z_old previous local solution
 * \param[in] dim block dimension
 * \param[in] omega relaxation parameter
 */
static inline void nsgs_relaxation(double* z_new, const double* z_old, int dim, double omega) {
  for (int d = 0; d < dim; ++d) {
    z_new[d] = omega * z_new[d] + (1.0 - omega) * z_old[d];
  }
}

/** Squared incremental error between two block vectors.
 *
 * \param[in] z_new new local solution
 * \param[in] z_old previous local solution
 * \param[in] dim block dimension
 * \return ||z_new - z_old||^2
 */
static inline double nsgs_light_error_squared(const double* z_new, const double* z_old,
                                              int dim) {
  double sum = 0.0;
  for (int d = 0; d < dim; ++d) {
    double diff = z_new[d] - z_old[d];
    sum += diff * diff;
  }
  return sum;
}

/** Squared Euclidean norm of a block vector.
 *
 * \param[in] z local solution vector
 * \param[in] dim block dimension
 * \return ||z||^2
 */
static inline double nsgs_squared_norm(const double* z, int dim) {
  double sum = 0.0;
  for (int d = 0; d < dim; ++d) {
    sum += z[d] * z[d];
  }
  return sum;
}

/** Determine convergence based on error and tolerance.
 *
 * \param[in] error current error
 * \param[in] tolerance convergence tolerance
 * \param[in] iter current iteration (unused, kept for API symmetry)
 * \param[in] options solver options (unused, kept for API symmetry)
 * \return 1 if not converged, 0 if converged
 */
static inline int nsgs_determine_convergence(double error, double tolerance, int iter,
                                             SolverOptions* options) {
  (void)iter;
  (void)options;
  return (error > tolerance) ? 1 : 0;
}

/** Check if too many blocks are frozen and reset if necessary.
 *
 * If nb_blocks - 1 or more blocks are frozen, reset all to unfreeze
 * at least one block for the next iteration.
 *
 * \param[in,out] freeze_blocks array of freeze counters
 * \param[in] nb_blocks number of blocks
 * \param[in] use_freezing whether freezing is enabled
 * \return number of frozen blocks
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

/** Default criteria for freezing a block (matches fc3d_nsgs).
 *
 * A block is frozen if either:
 * 1. Relative convergence: incremental_error <= tolerance^2 * 100^2 * ||z_local||^2
 * 2. Small solution: ||z_local||^2 <= ||z_global||^2 / (nb_blocks^2 * 1000)
 *
 * And we are past iteration 10.
 *
 * \param[in] incremental_error_2 squared incremental error for this block
 * \param[in] squared_norm_local squared norm of local solution
 * \param[in] tolerance global tolerance
 * \param[in] norm_z_global norm of global solution
 * \param[in] nb_blocks number of blocks
 * \param[in] iter current iteration
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

/** Update freeze counters for a block.
 *
 * \param[in,out] freeze_blocks array of freeze counters
 * \param[in] block block index
 * \param[in] freezing_iter number of iterations to freeze
 * \param[in] should_freeze whether this block should be frozen
 */
static inline void nsgs_update_freeze_counter(unsigned int* freeze_blocks, unsigned int block,
                                              int freezing_iter, int should_freeze) {
  if (!freeze_blocks) return;

  if (should_freeze) {
    freeze_blocks[block] = freezing_iter;
  }
}

/** Shuffle a block index array.
 *
 * \param[in,out] sblocks array of block indices to shuffle
 * \param[in] nb_blocks number of blocks
 * \param[in] shuffle_mode shuffle mode (0=none, 1=once, 2=each loop)
 * \param[in] iter current iteration
 */
static inline void nsgs_shuffle_blocks(unsigned int* sblocks, unsigned int nb_blocks,
                                       int shuffle_mode, int iter) {
  if (!sblocks) return;
  if (shuffle_mode == 0) return;

  /* Shuffle each loop, or first iteration only */
  if (shuffle_mode == 2 || (shuffle_mode == 1 && iter == 1)) {
    uint_shuffle(sblocks, nb_blocks);
  }
}

/** Count frozen blocks and return percentage.
 *
 * \param[in] freeze_blocks array of freeze counters
 * \param[in] nb_blocks number of blocks
 * \param[in] use_freezing whether freezing is enabled
 * \return percentage of frozen blocks (0-100)
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

/** Print iteration statistics (verbose output).
 *
 * \param[in] iter current iteration
 * \param[in] incremental_error incremental error value
 * \param[in] full_error full error value
 * \param[in] tolerance working tolerance
 * \param[in] user_tolerance user tolerance
 * \param[in] frozen_percent percentage of frozen blocks
 * \param[in] hasNotConverged convergence flag
 * \param[in] verbose verbosity level
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

/** Adapt working tolerance based on incremental and full error.
 *
 * When incremental error indicates convergence but full error does not,
 * adapt the working tolerance to try to achieve full convergence.
 *
 * \param[in] incremental_error current incremental error (squared norm of change)
 * \param[in] full_error current full error
 * \param[in] user_tolerance user-specified tolerance
 * \param[in,out] working_tolerance current working tolerance (may be adapted)
 * \param[in,out] localsolver_options internal solver options (may be tightened)
 * \param[in] iter current iteration
 * \param[in] verbose verbosity level
 * \return 0 if converged, 1 if not converged
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

#endif /* NSGS_UTILS_H */
