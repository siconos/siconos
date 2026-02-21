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

/*!\file NonSmoothGaussSeidel_options.h
 * \brief Options and constants for Non-Smooth Gauss-Seidel (NSGS) solvers
 *
 * This header contains all NSGS-specific options that are generic across
 * different problem formulations (friction contact, rolling friction, etc.)
 */

#ifndef NONSMOOTHGAUSSSEIDEL_OPTIONS_H
#define NONSMOOTHGAUSSSEIDEL_OPTIONS_H

/* ===========================================================================
 * NSGS Solver IDs
 * ===========================================================================
 * These are the main solver IDs for NSGS variants across different problem
 * formulations. They are kept in their respective formulation headers but
 * listed here for reference:
 *
 * Friction Contact 2D: SICONOS_FRICTION_2D_NSGS
 * Friction Contact 3D: SICONOS_FRICTION_3D_NSGS
 * Friction Contact 3D (velocity): SICONOS_FRICTION_3D_NSGSV
 * Global Friction 3D: SICONOS_GLOBAL_FRICTION_3D_NSGS
 * Global Friction 3D (WR): SICONOS_GLOBAL_FRICTION_3D_NSGS_WR
 * Rolling Friction 3D: SICONOS_ROLLING_FRICTION_3D_NSGS
 * Rolling Friction 2D: SICONOS_ROLLING_FRICTION_2D_NSGS
 */

/* ===========================================================================
 * NSGS Integer Parameters (iparam indices)
 * ===========================================================================
 * These indices are used in the iparam array of SolverOptions to configure
 * NSGS behavior. They are designed to coexist with other solver options.
 *
 * Note: Index 0 is reserved for SICONOS_IPARAM_MAX_ITER
 *       Index 1 is reserved for SICONOS_IPARAM_ITER_DONE
 */

/** NSGS parameter indices - iparam block starting at offset 10
 * This avoids conflict with general solver options (0-9)
 */
enum SICONOS_NSGS_IPARAM {
  /** Base offset for NSGS iparam to avoid conflicts */
  SICONOS_NSGS_IPARAM_OFFSET = 10,
  
  /** index in iparam to store the relaxation strategy */
  SICONOS_NSGS_RELAXATION = 10,
  
  /** index in iparam to store the shuffle strategy */
  SICONOS_NSGS_SHUFFLE = 11,
  
  /** index in iparam to store the shuffle seed */
  SICONOS_NSGS_SHUFFLE_SEED = 12,
  
  /** index in iparam to store the error evaluation method */
  SICONOS_NSGS_ERROR_EVALUATION_TYPE = 13,
  
  /** index in iparam to store the frequency of error evaluation */
  SICONOS_NSGS_ERROR_EVALUATION_FREQUENCY = 14,
  
  /** index in iparam to store the freezing contact strategy */
  SICONOS_NSGS_FREEZING_CONTACT = 15,
  
  /** index in iparam to store the filter local solution flag */
  SICONOS_NSGS_FILTER_LOCAL_SOLUTION = 16,
  
  /** index in iparam to store printing style (IPM-like) */
  SICONOS_NSGS_PRINTING_STYLE = 17,
  
  /** index in iparam for localsolver trivial solution flag */
  SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION = 18,
  
  /** index in iparam to store the current block/contact number */
  SICONOS_NSGS_CURRENT_BLOCK_NUMBER = 19
};

/* ===========================================================================
 * NSGS Double Parameters (dparam indices)
 * ===========================================================================
 * These indices are used in the dparam array of SolverOptions to configure
 * NSGS numerical parameters.
 *
 * Note: Index 0 is reserved for SICONOS_DPARAM_TOL
 *       Index 1 is reserved for SICONOS_DPARAM_RESIDU
 */

enum SICONOS_NSGS_DPARAM {
  /** Base offset for NSGS dparam to avoid conflicts */
  SICONOS_NSGS_DPARAM_OFFSET = 10,
  
  /** index in dparam to store the relaxation parameter omega */
  SICONOS_NSGS_RELAXATION_VALUE = 10,
  
  /** index in dparam to store the internal error ratio */
  SICONOS_NSGS_INTERNAL_ERROR_RATIO = 11
};

/* ===========================================================================
 * NSGS Error Evaluation Types
 * ===========================================================================
 * These enums define how the convergence error is evaluated during NSGS
 * iterations.
 */

enum SICONOS_NSGS_ERROR_EVALUATION {
  /** Evaluation of the error with the expensive full error function
   * This computes the true residual but is computationally expensive.
   */
  SICONOS_NSGS_ERROR_EVALUATION_FULL = 0,
  
  /** Evaluation of the error with the cheap incremental variation
   * Uses the norm of the difference between successive iterates.
   * Fast but may not reflect the true residual.
   */
  SICONOS_NSGS_ERROR_EVALUATION_LIGHT = 1,
  
  /** Evaluation with incremental variation but with full final check
   * Uses incremental error during iterations, but checks full error
   * when incremental error is below tolerance. Adapts tolerance if needed.
   */
  SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL = 2,
  
  /** Evaluation with adaptive frequency for full error computation
   * Starts with infrequent full error checks and increases frequency
   * as convergence is approached.
   */
  SICONOS_NSGS_ERROR_EVALUATION_ADAPTIVE = 3
};

/* ===========================================================================
 * NSGS Shuffle Strategies
 * ===========================================================================
 * These enums define how the order of block processing is shuffled.
 */

enum SICONOS_NSGS_SHUFFLE {
  /** No shuffling - process blocks in natural order */
  SICONOS_NSGS_SHUFFLE_FALSE = 0,
  
  /** Shuffle once at the beginning */
  SICONOS_NSGS_SHUFFLE_TRUE = 1,
  
  /** Reshuffle at the beginning of each iteration */
  SICONOS_NSGS_SHUFFLE_EACH_LOOP = 2
};

/* ===========================================================================
 * NSGS Relaxation Modes
 * ===========================================================================
 * These enums define whether relaxation is applied to the local solution.
 */

enum SICONOS_NSGS_RELAXATION {
  /** No relaxation */
  SICONOS_NSGS_RELAXATION_FALSE = 0,
  /** Apply relaxation with parameter omega */
  SICONOS_NSGS_RELAXATION_TRUE = 1
};

/* ===========================================================================
 * NSGS Local Solution Filtering
 * ===========================================================================
 * These enums define whether to filter (discard) local solutions with
 * suspicious residual values.
 */

enum SICONOS_NSGS_FILTER_LOCAL_SOLUTION {
  /** Do not filter local solutions */
  SICONOS_NSGS_FILTER_LOCAL_SOLUTION_FALSE = 0,
  /** Filter local solutions with NaN, Inf, or large residuals */
  SICONOS_NSGS_FILTER_LOCAL_SOLUTION_TRUE = 1
};

/* ===========================================================================
 * NSGS Printing Styles
 * ===========================================================================
 * These enums define the output printing format.
 */

enum SICONOS_NSGS_PRINTING_STYLE {
  /** Standard printing format */
  SICONOS_NSGS_PRINTING_STYLE_STANDARD = 0,
  /** IPM-like printing format */
  SICONOS_NSGS_PRINTING_STYLE_IPM = 1
};

/* ===========================================================================
 * NSGS Local Solver Options
 * ===========================================================================
 * These enums control the behavior of the internal local solver.
 */

enum SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION {
  /** Do not use trivial solution check */
  SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE = 0,
  /** Check for and use trivial solution when possible */
  SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE = 1
};

/* ===========================================================================
 * NSGS Internal Error Strategy
 * ===========================================================================
 * These enums define how the internal solver tolerance is adapted.
 */

enum SICONOS_NSGS_INTERNAL_ERROR_STRATEGY {
  /** Adaptive strategy based on current error */
  SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE = 0,
  /** Use a fixed given value */
  SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE = 1,
  /** Adaptive strategy based on number of contacts */
  SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT = 2
};

/* ===========================================================================
 * Backward Compatibility Aliases
 * ===========================================================================
 * These macros maintain compatibility with existing code that uses the old
 * Friction_cst.h naming convention.
 */

/* Old iparam names - mapped to new generic names */
#define SICONOS_FRICTION_3D_NSGS_RELAXATION           SICONOS_NSGS_RELAXATION
#define SICONOS_FRICTION_3D_NSGS_SHUFFLE              SICONOS_NSGS_SHUFFLE
#define SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED         SICONOS_NSGS_SHUFFLE_SEED
#define SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT     SICONOS_NSGS_FREEZING_CONTACT
#define SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION SICONOS_NSGS_FILTER_LOCAL_SOLUTION
#define SICONOS_FRICTION_3D_NSGS_PRINTING_LIKE_IPM    SICONOS_NSGS_PRINTING_STYLE

/* Old dparam names - mapped to new generic names */
#define SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE     SICONOS_NSGS_RELAXATION_VALUE

/* Old enum values - mapped to new generic names */
#define SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL                SICONOS_NSGS_ERROR_EVALUATION_FULL
#define SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT               SICONOS_NSGS_ERROR_EVALUATION_LIGHT
#define SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL
#define SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE            SICONOS_NSGS_ERROR_EVALUATION_ADAPTIVE

#define SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE        SICONOS_NSGS_SHUFFLE_FALSE
#define SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE         SICONOS_NSGS_SHUFFLE_TRUE
#define SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP SICONOS_NSGS_SHUFFLE_EACH_LOOP

#define SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE     SICONOS_NSGS_RELAXATION_FALSE
#define SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE      SICONOS_NSGS_RELAXATION_TRUE

#define SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE SICONOS_NSGS_FILTER_LOCAL_SOLUTION_FALSE
#define SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE  SICONOS_NSGS_FILTER_LOCAL_SOLUTION_TRUE

#define SICONOS_FRICTION_3D_NSGS_PRINTING_LIKE_IPM_FALSE SICONOS_NSGS_PRINTING_STYLE_STANDARD
#define SICONOS_FRICTION_3D_NSGS_PRINTING_LIKE_IPM_TRUE  SICONOS_NSGS_PRINTING_STYLE_IPM

#define SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE
#define SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE  SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE

/* Local solver iparam backward compatibility */
#define SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION

/* Error evaluation iparam backward compatibility */
#define SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION   SICONOS_NSGS_ERROR_EVALUATION_TYPE
#define SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY SICONOS_NSGS_ERROR_EVALUATION_FREQUENCY

/* Internal error strategy backward compatibility */
#define SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE       SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE
#define SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE    SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE
#define SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT

/* Old internal error ratio dparam */
#define SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO SICONOS_NSGS_INTERNAL_ERROR_RATIO

#endif /* NONSMOOTHGAUSSSEIDEL_OPTIONS_H */
