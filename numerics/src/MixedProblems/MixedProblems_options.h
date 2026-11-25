// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file MixedProblems_options.h
 * \brief Solver options and constants for Mixed Problems
 *
 * This header defines:
 * - Solver IDs (enum MIXED_SOLVER): All available solvers for mixed problems
 * - Integer parameter indices (iparam): Configuration options
 * - Double parameter indices (dparam): Numerical parameters
 *
 * Mixed problems combine heterogeneous model types (friction 2D/3D,
 * rolling friction 2D/3D, Von Mises plasticity, Drucker-Prager plasticity)
 * in a single unified problem structure with shared M matrix and q vector.
 */

#ifndef MIXEDPROBLEMS_OPTIONS_H
#define MIXEDPROBLEMS_OPTIONS_H

#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/* ===========================================================================
 * Solver IDs (enum MIXED_SOLVER)
 * ===========================================================================
 * Each solver is assigned a unique integer ID to avoid string-based errors.
 * IDs 700-799 are reserved for Mixed Problems solvers.
 */

/** \enum MIXED_SOLVER
 * \brief Unique identifiers for all Mixed Problem solvers
 */
enum MIXED_SOLVER {
  /* -----------------------------------------------------------------------
   * Mixed Problems solvers (IDs 700-799)
   * -----------------------------------------------------------------------
   * These solvers work on MixedProblems that combine heterogeneous
   * model types (friction, rolling friction, plasticity) in a single
   * unified problem structure.
   */

  /** Non-smooth Gauss-Seidel for MixedProblems */
  SICONOS_MIXED_NSGS = 700,

}; /* end of enum MIXED_SOLVER */

/* ===========================================================================
 * Integer Parameter Indices (iparam)
 * ===========================================================================
 * These enums define indices into the iparam array of SolverOptions.
 */

/** Indices for Mixed Problem solver integer parameters */
enum SICONOS_MIXED_IPARAM {
  /** Maximum number of iterations */
  SICONOS_MIXED_IPARAM_MAX_ITER = 0,

  /** Number of iterations done */
  SICONOS_MIXED_IPARAM_ITER_DONE = 1,

  /** Error evaluation method (see SICONOS_MIXED_ERROR_EVALUATION) */
  SICONOS_MIXED_IPARAM_ERROR_EVALUATION = 2,

  /** Frequency of error evaluation */
  SICONOS_MIXED_IPARAM_ERROR_EVALUATION_FREQUENCY = 3,

  /** Shuffle strategy (0=none, 1=once, 2=each iteration) */
  SICONOS_MIXED_IPARAM_SHUFFLE = 4,

  /** Current block number (for debugging) */
  SICONOS_MIXED_IPARAM_CURRENT_BLOCK = 5
};

/* ===========================================================================
 * Double Parameter Indices (dparam)
 * ===========================================================================
 */

/** Indices for Mixed Problem solver double parameters */
enum SICONOS_MIXED_DPARAM {
  /** Tolerance for convergence */
  SICONOS_MIXED_DPARAM_TOL = 0,

  /** Computed residual/error */
  SICONOS_MIXED_DPARAM_RESIDU = 1,

  /** Relaxation parameter (omega) */
  SICONOS_MIXED_DPARAM_RELAXATION = 2
};

/* ===========================================================================
 * Error Evaluation Types
 * ===========================================================================
 */

enum SICONOS_MIXED_ERROR_EVALUATION {
  /** Full error computation with expensive residual calculation */
  SICONOS_MIXED_ERROR_EVALUATION_FULL = 0,

  /** Light error using incremental variation */
  SICONOS_MIXED_ERROR_EVALUATION_LIGHT = 1,

  /** Light during iterations, full at the end */
  SICONOS_MIXED_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL = 2
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* MIXEDPROBLEMS_OPTIONS_H */
