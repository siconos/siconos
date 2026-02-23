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

/*!\file tolerance_manager.h
 * \brief Unified tolerance management for NSGS and related solvers
 *
 * This header provides a standardized way to handle:
 * - User-specified tolerance (original requirement)
 * - Working tolerance (adapted during iterations)
 * - Local solver tolerance (for internal solvers)
 * - Minimum tolerance bounds (numerical stability)
 */

#ifndef TOLERANCE_MANAGER_H
#define TOLERANCE_MANAGER_H

#include <float.h>
#include <math.h>
#include "SolverOptions.h"
#include "numerics_verbose.h"

/** Tolerance manager structure
 *
 * Encapsulates all tolerance-related state to avoid scattering
 * tolerance logic across solver implementations.
 */
typedef struct {
  /* User-specified parameters (read-only after init) */
  double user_tolerance;        /**< Original user-requested tolerance */
  double min_tolerance;         /**< Absolute floor (DBL_EPSILON * 1e-6) */
  
  /* Working state (may change during iterations) */
  double working_tolerance;     /**< Current adapted tolerance */
  double local_tol_original;    /**< Saved original local solver tolerance */
  double local_tol_current;     /**< Current (possibly adapted) local tolerance */
  
  /* Adaptation tracking */
  int adaptation_count;         /**< Number of times tolerance was adapted */
  double last_error_ratio;      /**< Last computed error ratio */
} ToleranceManager;

/** Initialize tolerance manager
 *
 * \param[in,out] tm Tolerance manager to initialize
 * \param[in] user_tol User-specified tolerance
 * \param[in] localsolver_options Local solver options (may be NULL)
 */
static inline void tolerance_manager_init(ToleranceManager* tm, 
                                          double user_tol,
                                          SolverOptions* localsolver_options) {
  tm->user_tolerance = user_tol;
  tm->working_tolerance = user_tol;
  tm->min_tolerance = DBL_EPSILON * 1e-6;
  tm->adaptation_count = 0;
  tm->last_error_ratio = 0.0;
  
  if (localsolver_options) {
    tm->local_tol_original = localsolver_options->dparam[SICONOS_DPARAM_TOL];
    tm->local_tol_current = tm->local_tol_original;
  } else {
    tm->local_tol_original = 0.0;
    tm->local_tol_current = 0.0;
  }
}

/** Get minimum allowed tolerance
 *
 * \param[in] tm Tolerance manager
 * \return Minimum tolerance (user_tol * 1e-6 or DBL_EPSILON * 1e-6)
 */
static inline double tolerance_manager_get_min(ToleranceManager* tm) {
  double user_min = tm->user_tolerance * 1e-6;
  return (user_min > tm->min_tolerance) ? user_min : tm->min_tolerance;
}

/** Adapt working tolerance based on error ratio
 *
 * When incremental error converges but full error doesn't, adapt the
 * working tolerance based on the error ratio.
 *
 * \param[in,out] tm Tolerance manager
 * \param[in] incr_error Incremental error (relative)
 * \param[in] full_error Full error (absolute)
 * \param[in] verbose Verbosity level
 * \return 0 if converged (full_error <= user_tolerance), 1 if not converged
 */
static inline int tolerance_manager_adapt_working(ToleranceManager* tm,
                                                  double incr_error,
                                                  double full_error,
                                                  int verbose) {
  /* Check if full error is within user tolerance */
  if (full_error <= tm->user_tolerance) {
    if (verbose > 0) {
      numerics_printf("ToleranceManager: CONVERGED - Full error (%.2e) <= user tolerance (%.2e)",
                      full_error, tm->user_tolerance);
    }
    return 0; /* Converged */
  }
  
  /* Not converged - compute adaptation */
  if (full_error > 0.0 && incr_error > 0.0) {
    double error_ratio = incr_error / full_error;
    double new_tolerance = error_ratio * tm->user_tolerance;
    double min_tol = tolerance_manager_get_min(tm);
    
    /* Cap at minimum tolerance */
    if (new_tolerance < min_tol) {
      if (verbose > 1) {
        numerics_printf("ToleranceManager: Capping tolerance at %.2e (computed: %.2e)",
                        min_tol, new_tolerance);
      }
      new_tolerance = min_tol;
    }
    
    /* Only adapt if tolerance would decrease */
    if (new_tolerance > 0.0 && new_tolerance < tm->working_tolerance) {
      if (verbose > 0) {
        numerics_printf("ToleranceManager: Adapting working tolerance: %.4e -> %.4e (ratio=%.4f)",
                        tm->working_tolerance, new_tolerance, error_ratio);
      }
      tm->working_tolerance = new_tolerance;
      tm->adaptation_count++;
      tm->last_error_ratio = error_ratio;
    }
  }
  
  return 1; /* Not converged */
}

/** Tighten local solver tolerance
 *
 * When incremental error is very small (NSGS not improving), tighten
 * the local solver tolerance instead of adapting working tolerance.
 *
 * \param[in,out] tm Tolerance manager
 * \param[in,out] localsolver_options Local solver options
 * \param[in] incr_error Current incremental error
 * \param[in] verbose Verbosity level
 * \return 1 (always continue iterating)
 */
static inline int tolerance_manager_tighten_local(ToleranceManager* tm,
                                                  SolverOptions* localsolver_options,
                                                  double incr_error,
                                                  int verbose) {
  if (!localsolver_options) {
    if (verbose > 0) {
      numerics_printf("ToleranceManager: Warning - No local solver options to tighten");
    }
    return 1;
  }
  
  double current_tol = tm->local_tol_current;
  double new_tol = fmax(current_tol / 100.0, tm->min_tolerance);
  
  localsolver_options->dparam[SICONOS_DPARAM_TOL] = new_tol;
  tm->local_tol_current = new_tol;
  
  if (verbose > 0) {
    numerics_printf("ToleranceManager: Incr error very small (%.2e), tightening local solver: %.2e -> %.2e",
                    incr_error, current_tol, new_tol);
  }
  
  return 1; /* Continue iterating */
}

/** Check if incremental error is "very small"
 *
 * \param[in] incr_error Incremental error
 * \return 1 if error is very small, 0 otherwise
 */
static inline int tolerance_manager_is_incr_very_small(double incr_error) {
  return (incr_error < 1e-15 || incr_error < DBL_EPSILON * 100.0) ? 1 : 0;
}

/** Handle convergence check with full error computation
 *
 * This is the main entry point for convergence checking when incremental
 * error has converged. It computes full error and either declares convergence
 * or adapts tolerance accordingly.
 *
 * \param[in,out] tm Tolerance manager
 * \param[in,out] localsolver_options Local solver options
 * \param[in] full_error Computed full error
 * \param[in] incr_error Incremental error
 * \param[in] verbose Verbosity level
 * \return 0 if converged, 1 if not converged
 */
static inline int tolerance_manager_check_convergence(ToleranceManager* tm,
                                                      SolverOptions* localsolver_options,
                                                      double full_error,
                                                      double incr_error,
                                                      int verbose) {
  /* First check: is full error within user tolerance? */
  if (full_error <= tm->user_tolerance) {
    if (verbose > 0) {
      numerics_printf("ToleranceManager: CONVERGED - Full (%.2e) <= user_tol (%.2e)",
                      full_error, tm->user_tolerance);
    }
    return 0; /* Converged */
  }
  
  /* Not converged - decide adaptation strategy */
  if (tolerance_manager_is_incr_very_small(incr_error)) {
    /* Incremental error negligible - tighten local solver */
    return tolerance_manager_tighten_local(tm, localsolver_options, incr_error, verbose);
  } else {
    /* Normal case - adapt working tolerance */
    return tolerance_manager_adapt_working(tm, incr_error, full_error, verbose);
  }
}

/** Restore original local solver tolerance
 *
 * Should be called at the end of the solver to restore the original
 * local solver tolerance.
 *
 * \param[in,out] tm Tolerance manager
 * \param[in,out] localsolver_options Local solver options
 */
static inline void tolerance_manager_restore_local(ToleranceManager* tm,
                                                   SolverOptions* localsolver_options) {
  if (localsolver_options) {
    localsolver_options->dparam[SICONOS_DPARAM_TOL] = tm->local_tol_original;
    tm->local_tol_current = tm->local_tol_original;
  }
}

/** Print tolerance manager state (for debugging)
 *
 * \param[in] tm Tolerance manager
 */
static inline void tolerance_manager_print(ToleranceManager* tm) {
  numerics_printf("ToleranceManager state:");
  numerics_printf("  User tolerance:     %.6e", tm->user_tolerance);
  numerics_printf("  Working tolerance:  %.6e", tm->working_tolerance);
  numerics_printf("  Minimum tolerance:  %.6e", tm->min_tolerance);
  numerics_printf("  Local tolerance:    %.6e (original: %.6e)",
                  tm->local_tol_current, tm->local_tol_original);
  numerics_printf("  Adaptation count:   %d", tm->adaptation_count);
}

#endif /* TOLERANCE_MANAGER_H */
