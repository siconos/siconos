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

/*!\file nsgs_local_problem.h
 * \brief Unified Local Problem Interface for NSGS
 *
 * This header provides a generic interface for local problem solvers used
 * within NSGS (Non-Smooth Gauss-Seidel) and other decomposition methods.
 *
 * Design Principles:
 * - Problem-type agnostic (works with FC3D, RFC3D, PLASTICITY_2D, etc.)
 * - Integrates with solver registration system
 * - Runtime local solver selection via registry
 */

#ifndef NSGS_LOCAL_PROBLEM_H
#define NSGS_LOCAL_PROBLEM_H

#include "SolverOptions.h"
#include "NonSmoothGaussSeidel_options.h"
#include "solver_registry.h"

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * \brief Local problem types
 */
typedef enum {
  NSGS_LP_FC3D = 0,   /**< 3D Friction Contact (dim 3) */
  NSGS_LP_FC2D,       /**< 2D Friction Contact (dim 2) */
  NSGS_LP_RFC3D,      /**< 3D Rolling Friction Contact (dim 5) */
  NSGS_LP_RFC2D,      /**< 2D Rolling Friction Contact (dim 3) */
  NSGS_LP_PLASTICITY_2D,       /**< 2D Mohr-Coulomb (dim 2) */
  NSGS_LP_MC3D        /**< 3D Mohr-Coulomb (dim 3) */
} NSGSLocalProblemType;

/**
 * \brief Opaque local problem handle
 */
typedef struct NSGSLocalProblem NSGSLocalProblem;

/**
 * \brief Local problem operations vtable
 */
typedef struct {
  /** Update local problem for given block */
  void (*update)(void* global, unsigned int block, const double* global_sol,
                 double* local_rhs, int dim, double* workspace);
  
  /** Extract local matrix for given block */
  void (*extract)(void* global, unsigned int block, double* local_mat, int dim);
  
  /** Solve local problem (optional, can use registry instead) */
  int (*solve)(void* local_data, const double* rhs, const double* mat,
               const double* mu, int dim, SolverOptions* opts, double* result);
  
  /** Problem type and dimension */
  NSGSLocalProblemType type;
  int dimension;
  
  /** Default local solver name */
  const char* default_solver;
} NSGSLocalProblemOps;

/* ============================================================================
 * Core API
 * ============================================================================ */

/**
 * \brief Create local problem handle
 */
NSGSLocalProblem* nsgs_local_problem_create(void* global_problem,
                                             unsigned int block_id,
                                             const NSGSLocalProblemOps* ops);

/**
 * \brief Free local problem
 */
void nsgs_local_problem_free(NSGSLocalProblem* local);

/**
 * \brief Update local problem with current global state
 */
void nsgs_local_problem_update(NSGSLocalProblem* local, const double* global_sol);

/**
 * \brief Solve local problem using registered solver
 */
int nsgs_local_problem_solve(NSGSLocalProblem* local, double* result,
                              SolverOptions* options);

/**
 * \brief Set local solver by name
 */
int nsgs_local_problem_set_solver(NSGSLocalProblem* local, const char* name);

/**
 * \brief Get problem dimension
 */
int nsgs_local_problem_get_dimension(const NSGSLocalProblem* local);

/**
 * \brief Get problem type
 */
NSGSLocalProblemType nsgs_local_problem_get_type(const NSGSLocalProblem* local);

/* ============================================================================
 * Problem-Specific Operations
 * ============================================================================ */

const NSGSLocalProblemOps* nsgs_fc3d_local_ops(void);
const NSGSLocalProblemOps* nsgs_fc2d_local_ops(void);
const NSGSLocalProblemOps* nsgs_rfc3d_local_ops(void);
const NSGSLocalProblemOps* nsgs_plasticity_2d_local_ops(void);

#if defined(__cplusplus)
}
#endif

#endif /* NSGS_LOCAL_PROBLEM_H */
