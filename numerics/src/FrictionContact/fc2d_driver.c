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

/*!\file fc2d_driver.c
 * \brief FC2D driver using the solver registration system
 */

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "fc2d_Solvers.h"
#include "fc3d_short_names.h"
#include "numerics_verbose.h"

/* Registration-based headers */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"

#include <stdio.h>
#include <stdlib.h>

/* ===========================================================================
 * Registration-Based Driver
 * =========================================================================== */

int fc2d_driver(FrictionContactProblem* problem, double* reaction,
                double* velocity, SolverOptions* options) {
  /* Input validation using standardized macros */
  CHECK_NULL(problem);
  CHECK_NULL(reaction);
  CHECK_NULL(velocity);
  CHECK_OPTIONS(options);
  CHECK_MATRIX(problem->M);
  CHECK_NULL(problem->q);
  CHECK_NULL(problem->mu);
  CHECK_ARG(problem->numberOfContacts > 0, "Number of contacts must be positive");

  /* Check dimension */
  CHECK_DIMENSION(problem->dimension, 2);

  /* Initialize output */
  SET_SOLVER_ITER_DONE(options, 0);
  SET_SOLVER_RESIDUAL(options, 0.0);

  /* Handle sparse matrices - convert to dense for non-NSGS solvers */
  NumericsMatrix* M_original = NULL;
  if (problem->M->storageType != NM_DENSE && options->solverId != FC2D_NSGS) {
    numerics_printf_verbose(1, "fc2d_driver: converting sparse matrix to dense for solver %d",
                            options->solverId);
    M_original = problem->M;
    problem->M = NM_create(NM_DENSE, M_original->size0, M_original->size1);
    NM_to_dense(M_original, problem->M);
  }

  /* Lookup solver in registry */
  const SolverEntry* solver = solver_registry_lookup(options->solverId);
  if (!solver) {
    numerics_printf("fc2d_driver: solver ID %d not found in registry", options->solverId);
    /* Restore original matrix if converted */
    if (M_original) {
      NM_free(problem->M);
      problem->M = M_original;
    }
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  numerics_printf_verbose(1, "fc2d_driver: using solver '%s' (%s)",
                          solver->name, solver->description);

  /* Validate solver is appropriate for this problem type */
  if (solver->is_local_solver) {
    numerics_printf("fc2d_driver: solver '%s' is a local solver, cannot be used as main solver",
                    solver->name);
    /* Restore original matrix if converted */
    if (M_original) {
      NM_free(problem->M);
      problem->M = M_original;
    }
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  /* Check solve function exists */
  if (!solver->solve) {
    numerics_printf("fc2d_driver: solver '%s' has no solve function", solver->name);
    /* Restore original matrix if converted */
    if (M_original) {
      NM_free(problem->M);
      problem->M = M_original;
    }
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  /* Initialize solver if init function provided */
  if (solver->init) {
    int init_status = solver->init(problem, options);
    if (init_status != NUMERICS_OK) {
      numerics_printf("fc2d_driver: solver initialization failed with error %d", init_status);
      /* Restore original matrix if converted */
      if (M_original) {
        NM_free(problem->M);
        problem->M = M_original;
      }
      return init_status;
    }
  }

  /* Call the solver */
  numerics_printf_verbose(1, "fc2d_driver: calling solver...");
  int solve_status = solver->solve(problem, reaction, velocity, options);

  /* Cleanup if needed */
  if (solver->free) {
    solver->free(problem, options);
  }

  /* Restore original matrix if converted */
  if (M_original) {
    NM_free(problem->M);
    problem->M = M_original;
  }

  /* Log result */
  if (solve_status == NUMERICS_OK) {
    numerics_printf_verbose(1, "fc2d_driver: solver converged successfully");
  } else {
    numerics_printf_verbose(1, "fc2d_driver: solver returned status %d", solve_status);
  }

  return solve_status;
}

/* ===========================================================================
 * Convenience Functions
 * =========================================================================== */

SolverOptions* fc2d_solver_options_create(solver_id_t solver_id) {
  const SolverEntry* solver = solver_registry_lookup(solver_id);
  
  if (!solver) {
    fprintf(stderr, "fc2d_solver_options_create: solver ID %d not registered\n", solver_id);
    return NULL;
  }
  
  if (solver->is_local_solver) {
    fprintf(stderr, "fc2d_solver_options_create: solver '%s' is a local solver\n", solver->name);
    return NULL;
  }
  
  SolverOptions* options = solver_options_create(solver_id);
  if (!options) {
    fprintf(stderr, "fc2d_solver_options_create: failed to create options\n");
    return NULL;
  }
  
  if (solver->init) {
    int init_status = solver->init(NULL, options);
    if (init_status != NUMERICS_OK) {
      fprintf(stderr, "fc2d_solver_options_create: init failed with status %d\n", init_status);
      solver_options_delete(options);
      return NULL;
    }
  }
  
  return options;
}

void fc2d_list_available_solvers(void) {
  printf("\nAvailable FC2D Solvers:\n");
  printf("%-8s %-20s %-12s %-12s\n", "ID", "Name", "Max Iter", "Tolerance");
  printf("%-8s %-20s %-12s %-12s\n", "--------", "--------------------", "------------", "------------");
  
  size_t count;
  const SolverEntry** solvers = solver_registry_get_all(&count);
  
  for (size_t i = 0; i < count; i++) {
    const SolverEntry* s = solvers[i];
    if (s->id >= 400 && s->id < 500) {
      printf("%-8d %-20s %-12d %-12.2e\n",
             s->id, s->name, s->default_max_iter, s->default_tol);
    }
  }
  printf("\n");
}

void fc2d_print_solver_info(solver_id_t solver_id) {
  const SolverEntry* solver = solver_registry_lookup(solver_id);
  
  if (!solver) {
    printf("Solver ID %d not found in registry.\n", solver_id);
    return;
  }
  
  printf("\nSolver: %s (ID: %d)\n", solver->name, solver->id);
  printf("Description: %s\n", solver->description);
  printf("Max iterations: %d\n", solver->default_max_iter);
  printf("Tolerance: %.2e\n", solver->default_tol);
  printf("\n");
}
