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

/*!\file naming_conventions.h
 * \brief Standardized naming conventions for numerics solvers
 *
 * This header defines consistent naming patterns and accessor macros
 * for variables across all numerics solvers.
 *
 * == DOMAIN-SPECIFIC VARIABLE NAMING ==
 *
 * FRICTION CONTACT (FrictionContact/, RollingFrictionContact/):
 *   Use physics-based names that match the problem domain:
 *   - 'reaction'    - Contact reaction vector (global)
 *   - 'velocity'    - Relative velocity vector (global)
 *   - 'r_local'     - Local reaction for a single contact
 *   - 'v_local'     - Local velocity for a single contact
 *
 * PLASTICITY (Plasticity/):
 *   Use mechanics-based names:
 *   - 'stress'      - Stress vector (global)
 *   - 'strainrate'  - Strain rate vector (global)
 *   - 's_local'     - Local stress for a single cone
 *   - 'e_local'     - Local strain rate for a single cone
 *
 * GENERIC IMPLEMENTATIONS (NonSmoothSolvers/, GenericMechanical/):
 *   Use mathematical notation:
 *   - 'x'           - Primal variable (e.g., reaction, stress)
 *   - 'z'           - Dual variable (e.g., velocity, strainrate)
 *   - 'x_local'     - Local primal variable
 *   - 'z_local'     - Local dual variable
 *
 * == STANDARDIZED TERMINOLOGY ==
 *
 * SUBSYSTEM IDENTIFIERS:
 * - 'block'         - Generic term for contact/cone/subsystem (preferred)
 * - 'contact'       - 3D friction contact (FC3D)
 * - 'cone'          - Mohr-Coulomb cone (MC2D) or friction cone
 * - 'nc' or 'n_blocks' - Number of subsystems
 * - 'block_id'      - Index of current subsystem (preferred over 'contact', 'cone')
 *
 * SOLVER OPTIONS - USE MACROS (not direct array access):
 * - SOLVER_TOL(options)        -> options->dparam[SICONOS_DPARAM_TOL]
 * - SOLVER_RESIDUAL(options)   -> options->dparam[SICONOS_DPARAM_RESIDU]
 * - SOLVER_MAX_ITER(options)   -> options->iparam[SICONOS_IPARAM_MAX_ITER]
 * - SOLVER_ITER_DONE(options)  -> options->iparam[SICONOS_IPARAM_ITER_DONE]
 *
 * SOLVER OPTIONS:
 * - 'options'       - Main solver options
 * - 'local_opts'    - Local/internal solver options (preferred over 'localsolver_options')
 * - 'user_tol'      - User-specified tolerance
 * - 'work_tol'      - Working tolerance (may be adapted)
 * - 'local_tol'     - Local solver tolerance
 *
 * ERRORS:
 * - 'err_full'      - Full/residual error
 * - 'err_incr'      - Incremental/light error
 * - 'err_rel'       - Relative error
 * - 'err_local'     - Local subsystem error
 *
 * ITERATION COUNTERS:
 * - 'iter'          - Current iteration number
 * - 'itermax'       - Maximum iterations
 * - 'iter_done'     - Iterations completed
 * - 'newton_iter'   - Newton iteration counter
 * - 'ls_iter'       - Line search iteration counter
 *
 * PROBLEM DATA:
 * - 'problem'       - Global problem structure
 * - 'localpb' or 'pb_local' - Local problem structure
 * - 'M'             - Global matrix
 * - 'M_local'       - Local subsystem matrix
 * - 'q'             - Global right-hand side
 * - 'q_local'       - Local right-hand side
 */

#ifndef NAMING_CONVENTIONS_H
#define NAMING_CONVENTIONS_H

/* ============================================================================
 * ACCESSOR MACROS FOR SOLVER OPTIONS
 * ============================================================================ */

/** Get user tolerance from solver options */
#define SOLVER_TOL(options) ((options)->dparam[SICONOS_DPARAM_TOL])

/** Get local solver tolerance */
#define LOCAL_SOLVER_TOL(local_opts) ((local_opts)->dparam[SICONOS_DPARAM_TOL])

/** Set local solver tolerance */
#define SET_LOCAL_SOLVER_TOL(local_opts, value) \
  ((local_opts)->dparam[SICONOS_DPARAM_TOL] = (value))

/** Get residual from solver options */
#define SOLVER_RESIDUAL(options) ((options)->dparam[SICONOS_DPARAM_RESIDU])

/** Set residual in solver options */
#define SET_SOLVER_RESIDUAL(options, value) \
  ((options)->dparam[SICONOS_DPARAM_RESIDU] = (value))

/** Get max iterations from solver options */
#define SOLVER_MAX_ITER(options) ((options)->iparam[SICONOS_IPARAM_MAX_ITER])

/** Get iterations done from solver options */
#define SOLVER_ITER_DONE(options) ((options)->iparam[SICONOS_IPARAM_ITER_DONE])

/** Set iterations done in solver options */
#define SET_SOLVER_ITER_DONE(options, value) \
  ((options)->iparam[SICONOS_IPARAM_ITER_DONE] = (value))

/** Get current block/contact number from solver options */
#define SOLVER_CURRENT_BLOCK(options) \
  ((options)->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER])

/** Set current block/contact number in solver options */
#define SET_SOLVER_CURRENT_BLOCK(options, block_id) \
  ((options)->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] = (block_id))

/* ============================================================================
 * COMMON VARIABLE SHORTCUTS
 * ============================================================================ */

/** Block size for 2D problems */
#define BLOCK_SIZE_2D 2

/** Block size for 3D problems */
#define BLOCK_SIZE_3D 3

/** Block size for 5D rolling friction problems */
#define BLOCK_SIZE_5D 5

/** Compute global index from block index and local dimension */
#define GLOBAL_INDEX(block_id, local_dim, local_idx) ((block_id) * (local_dim) + (local_idx))

/** Compute block start index */
#define BLOCK_START(block_id, local_dim) ((block_id) * (local_dim))

/* ============================================================================
 * CONSISTENT TYPE ALIASES
 * ============================================================================ */

/** Type for block/contact indices - use unsigned for array indices */
typedef unsigned int block_id_t;

/** Type for iteration counters */
typedef int iter_t;

/** Type for error values */
typedef double numerical_error_t;

/** Type for tolerance values */
typedef double tolerance_t;

/* ============================================================================
 * DEPRECATED NAMES - Do not use in new code
 * ============================================================================ */

/* These are kept for backward compatibility but should not be used in new code:
 *
 * INCONSISTENT            PREFERRED
 * ---------               ---------
 * contact (as index)  ->  block_id
 * cone (as index)     ->  block_id
 * R (reaction)        ->  reaction or z
 * var_z               ->  z or reaction
 * var_x               ->  x or velocity
 * localsolver_options ->  local_opts
 * dparam[0]           ->  dparam[SICONOS_DPARAM_TOL]
 * iparam[0]           ->  iparam[SICONOS_IPARAM_MAX_ITER]
 */

#endif /* NAMING_CONVENTIONS_H */
