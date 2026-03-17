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

#include "nsgs_local_problem.h"
#include "numerics_verbose.h"
#include "numerics_errors.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Internal Structure
 * ============================================================================ */

struct NSGSLocalProblem {
  void* global_problem;
  unsigned int block_id;
  NSGSLocalProblemOps ops;
  
  /* Buffers */
  double* local_rhs;
  double* local_mat;
  double* workspace;
  
  /* Solver */
  const SolverEntry* solver;
  SolverOptions* solver_opts;
};

/* ============================================================================
 * Core Implementation
 * ============================================================================ */

NSGSLocalProblem* nsgs_local_problem_create(void* global_problem,
                                             unsigned int block_id,
                                             const NSGSLocalProblemOps* ops) {
  if (!global_problem || !ops) {
    numerics_error("nsgs_local_problem_create", "NULL arguments");
    return NULL;
  }
  
  NSGSLocalProblem* local = (NSGSLocalProblem*)calloc(1, sizeof(NSGSLocalProblem));
  if (!local) return NULL;
  
  local->global_problem = global_problem;
  local->block_id = block_id;
  local->ops = *ops;
  
  int dim = ops->dimension;
  local->local_rhs = (double*)calloc(dim, sizeof(double));
  local->local_mat = (double*)calloc(dim * dim, sizeof(double));
  local->workspace = (double*)calloc(dim * dim, sizeof(double));
  
  if (!local->local_rhs || !local->local_mat || !local->workspace) {
    nsgs_local_problem_free(local);
    return NULL;
  }
  
  return local;
}

void nsgs_local_problem_free(NSGSLocalProblem* local) {
  if (!local) return;
  free(local->local_rhs);
  free(local->local_mat);
  free(local->workspace);
  if (local->solver_opts) {
    solver_options_delete(local->solver_opts);
  }
  free(local);
}

void nsgs_local_problem_update(NSGSLocalProblem* local, const double* global_sol) {
  if (!local || !global_sol) return;
  
  if (local->ops.update) {
    local->ops.update(local->global_problem, local->block_id, global_sol,
                      local->local_rhs, local->ops.dimension, local->workspace);
  }
  
  if (local->ops.extract) {
    local->ops.extract(local->global_problem, local->block_id,
                       local->local_mat, local->ops.dimension);
  }
}

int nsgs_local_problem_solve(NSGSLocalProblem* local, double* result,
                              SolverOptions* options) {
  if (!local || !result) {
    return NUMERICS_ERR_NULL_POINTER;
  }
  
  /* Use provided options or local solver options */
  SolverOptions* opts = options ? options : local->solver_opts;
  
  /* If custom solve function provided, use it */
  if (local->ops.solve) {
    return local->ops.solve(local, local->local_rhs, local->local_mat, NULL,
                            local->ops.dimension, opts, result);
  }
  
  /* Otherwise use registered solver */
  if (!local->solver) {
    numerics_error("nsgs_local_problem_solve", "No local solver configured");
    return NUMERICS_ERR_INVALID_SOLVER;
  }
  
  if (!local->solver->solve) {
    return NUMERICS_ERR_INVALID_SOLVER;
  }
  
  /* Call registered solver - local problem passed as problem data */
  return local->solver->solve(local->global_problem, result, NULL, opts);
}

int nsgs_local_problem_set_solver(NSGSLocalProblem* local, const char* name) {
  if (!local || !name) {
    return NUMERICS_ERR_NULL_POINTER;
  }
  
  const SolverEntry* solver = solver_registry_lookup_by_name(name);
  CHECK_COND(solver != NULL, NUMERICS_ERR_INVALID_SOLVER, "Solver not found");
  
  local->solver = solver;
  
  /* Create/update solver options */
  if (local->solver_opts) {
    solver_options_delete(local->solver_opts);
  }
  local->solver_opts = solver_options_create(solver->id);
  if (local->solver_opts) {
    local->solver_opts->iparam[SICONOS_IPARAM_MAX_ITER] = solver->default_max_iter;
    local->solver_opts->dparam[SICONOS_DPARAM_TOL] = solver->default_tol;
  }
  
  return NUMERICS_OK;
}

int nsgs_local_problem_get_dimension(const NSGSLocalProblem* local) {
  return local ? local->ops.dimension : -1;
}

NSGSLocalProblemType nsgs_local_problem_get_type(const NSGSLocalProblem* local) {
  return local ? local->ops.type : NSGS_LP_FC3D;
}

/* ============================================================================
 * Problem-Specific Operations (Stubs - implement in problem files)
 * ============================================================================ */

static void default_update(void* global, unsigned int block, const double* global_sol,
                           double* local_rhs, int dim, double* workspace) {
  (void)global; (void)block; (void)global_sol; (void)local_rhs; (void)dim; (void)workspace;
}

static void default_extract(void* global, unsigned int block, double* local_mat, int dim) {
  (void)global; (void)block; (void)local_mat; (void)dim;
}

static const NSGSLocalProblemOps fc3d_ops = {
  .update = default_update,
  .extract = default_extract,
  .solve = NULL,
  .type = NSGS_LP_FC3D,
  .dimension = 3,
  .default_solver = "PROJECTION"
};

const NSGSLocalProblemOps* nsgs_fc3d_local_ops(void) {
  return &fc3d_ops;
}

static const NSGSLocalProblemOps fc2d_ops = {
  .update = default_update,
  .extract = default_extract,
  .solve = NULL,
  .type = NSGS_LP_FC2D,
  .dimension = 2,
  .default_solver = "PROJECTION"
};

const NSGSLocalProblemOps* nsgs_fc2d_local_ops(void) {
  return &fc2d_ops;
}

static const NSGSLocalProblemOps rfc3d_ops = {
  .update = default_update,
  .extract = default_extract,
  .solve = NULL,
  .type = NSGS_LP_RFC3D,
  .dimension = 5,
  .default_solver = "PROJECTION"
};

const NSGSLocalProblemOps* nsgs_rfc3d_local_ops(void) {
  return &rfc3d_ops;
}

static const NSGSLocalProblemOps plasticity_2d_ops = {
  .update = default_update,
  .extract = default_extract,
  .solve = NULL,
  .type = NSGS_LP_PLASTICITY_2D,
  .dimension = 2,
  .default_solver = "PROJECTION"
};

const NSGSLocalProblemOps* nsgs_plasticity_2d_local_ops(void) {
  return &plasticity_2d_ops;
}
