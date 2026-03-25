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

/*!\file error_computation.h
 * \brief Standardized error computation for iterative solvers
 *
 * This header provides consistent error computation functions for:
 * - Incremental/light error (||z_new - z_old||)
 * - Full/residual error (problem-dependent)
 * - Relative error normalization
 */

#ifndef ERROR_COMPUTATION_H
#define ERROR_COMPUTATION_H

#include <float.h>
#include <math.h>
#include "SiconosBlas.h"

/** Compute incremental (light) error
 *
 * Computes sqrt(sum of squared differences between new and old values).
 * This is efficient for checking NSGS convergence between iterations.
 *
 * \param[in] z_new New solution vector
 * \param[in] z_old Old solution vector
 * \param[in] dim Dimension of vectors
 * \param[in] normalize Whether to normalize by norm of z_new
 * \return Incremental error (relative if normalize=1, absolute otherwise)
 */
static inline double error_compute_incremental(const double* z_new,
                                               const double* z_old,
                                               unsigned int dim,
                                               int normalize) {
  double error_sum = 0.0;
  double norm_sum = 0.0;
  
  for (unsigned int i = 0; i < dim; ++i) {
    double diff = z_new[i] - z_old[i];
    error_sum += diff * diff;
    if (normalize) {
      norm_sum += z_new[i] * z_new[i];
    }
  }
  
  double error = sqrt(error_sum);
  
  if (normalize && norm_sum > DBL_EPSILON) {
    error /= sqrt(norm_sum);
  }
  
  return error;
}

/** Compute incremental error using BLAS (for large vectors)
 *
 * \param[in] z_new New solution vector
 * \param[in] z_old Old solution vector
 * \param[in] dim Dimension of vectors
 * \param[in] normalize Whether to normalize by norm of z_new
 * \return Incremental error
 */
static inline double error_compute_incremental_blas(const double* z_new,
                                                    const double* z_old,
                                                    unsigned int dim,
                                                    int normalize) {
  /* Use dnrm2 for diff = z_new - z_old */
  double* diff = (double*)malloc(dim * sizeof(double));
  if (!diff) return -1.0; /* Error */
  
  for (unsigned int i = 0; i < dim; ++i) {
    diff[i] = z_new[i] - z_old[i];
  }
  
  double error = cblas_dnrm2((int32_t)dim, diff, 1);
  
  if (normalize) {
    double norm = cblas_dnrm2((int32_t)dim, z_new, 1);
    if (norm > DBL_EPSILON) {
      error /= norm;
    }
  }
  
  free(diff);
  return error;
}

/** Compute relative error
 *
 * \param[in] error Absolute error
 * \param[in] reference Reference value (norm of solution)
 * \return Relative error
 */
static inline double error_make_relative(double error, double reference) {
  if (reference > DBL_EPSILON) {
    return error / reference;
  }
  return error;
}

/** Check if error is converged within tolerance
 *
 * \param[in] error Current error
 * \param[in] tolerance Convergence tolerance
 * \return 1 if converged, 0 otherwise
 */
static inline int error_is_converged(double error, double tolerance) {
  return (error <= tolerance) ? 1 : 0;
}

/** Compute squared norm of difference (for use in incremental error accumulation)
 *
 * \param[in] z_new New solution vector
 * \param[in] z_old Old solution vector
 * \param[in] dim Dimension
 * \return Sum of squared differences
 */
static inline double error_incr_squared_sum(const double* z_new,
                                            const double* z_old,
                                            unsigned int dim) {
  double sum = 0.0;
  for (unsigned int i = 0; i < dim; ++i) {
    double diff = z_new[i] - z_old[i];
    sum += diff * diff;
  }
  return sum;
}

/** Finalize incremental error computation
 *
 * Converts sum of squared differences to relative error.
 *
 * \param[in] squared_sum Sum of squared differences
 * \param[in] dim Total dimension
 * \param[in] z_full Full solution vector (for normalization)
 * \return Relative incremental error
 */
static inline double error_incr_finalize(double squared_sum,
                                         unsigned int dim,
                                         const double* z_full) {
  double error = sqrt(squared_sum);
  double norm = cblas_dnrm2((int32_t)dim, z_full, 1);
  
  if (norm > DBL_EPSILON) {
    error /= norm;
  }
  
  return error;
}

/** Error statistics structure for tracking convergence history */
typedef struct {
  double initial_error;    /**< Error at first iteration */
  double previous_error;   /**< Error at previous iteration */
  double current_error;    /**< Error at current iteration */
  double best_error;       /**< Best error seen so far */
  int stagnation_count;    /**< Number of iterations with small improvement */
} ErrorStats;

/** Initialize error statistics
 *
 * \param[in,out] stats Error statistics structure
 */
static inline void error_stats_init(ErrorStats* stats) {
  stats->initial_error = -1.0;
  stats->previous_error = -1.0;
  stats->current_error = -1.0;
  stats->best_error = 1e300;
  stats->stagnation_count = 0;
}

/** Update error statistics
 *
 * \param[in,out] stats Error statistics structure
 * \param[in] new_error New error value
 */
static inline void error_stats_update(ErrorStats* stats, double new_error) {
  stats->previous_error = stats->current_error;
  stats->current_error = new_error;
  
  if (stats->initial_error < 0.0) {
    stats->initial_error = new_error;
  }
  
  if (new_error < stats->best_error) {
    stats->best_error = new_error;
    stats->stagnation_count = 0;
  } else {
    /* Check for stagnation (less than 1% improvement) */
    if (stats->previous_error > 0.0) {
      double improvement = (stats->previous_error - new_error) / stats->previous_error;
      if (improvement < 0.01) {
        stats->stagnation_count++;
      }
    }
  }
}

/** Check if solver is stagnating
 *
 * \param[in] stats Error statistics
 * \param[in] threshold Number of iterations to consider as stagnation
 * \return 1 if stagnating, 0 otherwise
 */
static inline int error_stats_is_stagnating(ErrorStats* stats, int threshold) {
  return (stats->stagnation_count >= threshold) ? 1 : 0;
}

/** Compute error reduction rate
 *
 * \param[in] stats Error statistics
 * \return Reduction rate (current / initial), -1.0 if not available
 */
static inline double error_stats_reduction_rate(ErrorStats* stats) {
  if (stats->initial_error > 0.0 && stats->current_error >= 0.0) {
    return stats->current_error / stats->initial_error;
  }
  return -1.0;
}

#endif /* ERROR_COMPUTATION_H */
