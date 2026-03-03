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

/*!\file fc3d_driver.c
 * \brief FC3D driver using the solver registration system
 *
 * This driver uses the solver registration system for dynamic dispatch.
 * No giant switch statement needed - just lookup and dispatch!
 *
 * COMPARISON:
 * - Old driver: ~200 lines of switch statements
 * - New driver: ~80 lines of clean, maintainable code
 *
 * BENEFITS:
 * - No hardcoded solver IDs in the driver
 * - Automatic support for new solvers (just register them)
 * - Consistent error handling
 * - Runtime solver introspection
 */

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "numerics_verbose.h"

/* Registration-based headers */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* Forward declarations for convenience functions */
void fc3d_list_available_solvers(void);
void fc3d_print_solver_info(solver_id_t solver_id);

/* ===========================================================================
 * Trivial Case Check
 * =========================================================================== */

int fc3d_checkTrivialCase(FrictionContactProblem* problem, double* velocity, double* reaction,
                          SolverOptions* options) {
  (void)options;
  assert(problem);
  assert(problem->q);
  
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  int n = 3 * nc;
  
  for (int i = 0; i < nc; i++) {
    if (q[3 * i] < -DBL_EPSILON) return NUMERICS_ERR_INFEASIBLE;
  }
  for (int i = 0; i < n; ++i) {
    velocity[i] = q[i];
    reaction[i] = 0.;
  }
  
  numerics_printf("fc3d: trivial solution (take-off), reaction = 0, velocity = q.");
  return NUMERICS_OK;
}

/* ===========================================================================
 * Main Driver - Registration-Based
 * =========================================================================== */

int fc3d_driver(FrictionContactProblem* problem, double* reaction, double* velocity,
                SolverOptions* options) {
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
  CHECK_DIMENSION(problem->dimension, 3);

  /* Initialize output */
  SET_SOLVER_ITER_DONE(options, 0);
  SET_SOLVER_RESIDUAL(options, 0.0);

  /* Check for trivial case */
  int trivial_status = fc3d_checkTrivialCase(problem, velocity, reaction, options);
  if (trivial_status == NUMERICS_OK) {
    return NUMERICS_OK;
  }

  /* Lookup solver in registry */
  const SolverEntry* solver = solver_registry_lookup(options->solverId);
  CHECK_COND(solver != NULL, NUMERICS_ERR_INVALID_SOLVER, 
             "Solver ID not found in registry");

  numerics_printf_verbose(1, "fc3d_driver: using solver '%s' (%s)",
                          solver->name, solver->description);

  /* /\* Validate solver is appropriate for this problem type *\/ */
  /* CHECK_COND(!solver->is_local_solver, NUMERICS_ERR_INVALID_SOLVER, */
  /*            "Local solver cannot be used as main solver"); */

  /* Check solve function exists */
  CHECK_COND(solver->solve != NULL, NUMERICS_ERR_INVALID_SOLVER,
             "Solver has no solve function");

  /* Initialize solver if init function provided */
  if (solver->init) {
    int init_status = solver->init(problem, options);
    if (init_status != NUMERICS_OK) {
      fprintf(stderr, "[ERROR] fc3d_driver: solver initialization failed: %s\n",
              numerics_error_string(init_status));
      return init_status;
    }
  }

  /* Call the solver */
  numerics_printf_verbose(1, "fc3d_driver: calling solver...");
  int solve_status = solver->solve(problem, reaction, velocity, options);

  /* Cleanup if needed */
  if (solver->free) {
    solver->free(problem, options);
  }

  /* Log result */
  if (solve_status == NUMERICS_OK) {
    numerics_printf_verbose(1, "fc3d_driver: solver converged successfully");
  } else {
    fprintf(stderr, "[WARN] fc3d_driver: solver returned status %d (%s)\n",
            solve_status, numerics_error_string(solve_status));
  }

  return solve_status;
}

/* ===========================================================================
 * Convenience Function: Create Options with Validation
 * =========================================================================== */

SolverOptions* fc3d_solver_options_create(solver_id_t solver_id) {
  /* Lookup solver first to validate it exists and is appropriate */
  const SolverEntry* solver = solver_registry_lookup(solver_id);
  
  if (!solver) {
    fprintf(stderr, "[ERROR] fc3d_solver_options_create: solver ID %d not registered\n", 
            solver_id);
    fprintf(stderr, "[INFO] Available FC3D solvers:\n");
    fc3d_list_available_solvers();
    return NULL;
  }
  
  if (solver->is_local_solver) {
    fprintf(stderr, "[ERROR] fc3d_solver_options_create: solver '%s' is a local solver, "
                    "use it within NSGS instead\n", solver->name);
    return NULL;
  }
  
  /* Create options using registered defaults */
  SolverOptions* options = solver_options_create(solver_id);
  
  if (!options) {
    fprintf(stderr, "[ERROR] fc3d_solver_options_create: failed to create options\n");
    return NULL;
  }
  
  /* Initialize with solver defaults if available */
  if (solver->init) {
    int init_status = solver->init(NULL, options);
    if (init_status != NUMERICS_OK) {
      fprintf(stderr, "fc3d_solver_options_create: init failed with status %d\n", 
              init_status);
      solver_options_delete(options);
      return NULL;
    }
  }
  
  return options;
}

/* ===========================================================================
 * Convenience Function: List Available Solvers
 * =========================================================================== */

void fc3d_list_available_solvers(void) {
  printf("\n");
  printf("+=============================================================================+\n");
  printf("| Available FC3D Solvers                                                       |\n");
  printf("+=============================================================================+\n");
  printf("| %-20s | %-8s | %-12s | %-12s | %-5s |\n", 
         "Name", "ID", "Max Iter", "Tolerance", "Type");
  printf("+----------------------+----------+--------------+--------------+-------+\n");
  
  size_t count;
  const SolverEntry** solvers = solver_registry_get_all(&count);
  
  int fc3d_count = 0;
  for (size_t i = 0; i < count; i++) {
    const SolverEntry* s = solvers[i];
    /* Show FC3D solvers (ID range 500-599 for 3D friction contact) */
    if ((s->id >= 500 && s->id < 600) || 
        (s->id >= 5000 && s->id < 6000)) {  /* Global formulation range */
      printf("| %-20s | %-8d | %-12d | %-12.2e | %-5s |\n",
             s->name, s->id, s->default_max_iter, s->default_tol,
             s->is_local_solver ? "local" : "main");
      fc3d_count++;
    }
  }
  
  printf("+----------------------+----------+--------------+--------------+-------+\n");
  printf("| Total: %d solver(s) available                                               |\n", 
         fc3d_count);
  printf("+=============================================================================+\n");
  printf("\n");
  printf("Usage: SolverOptions* options = fc3d_solver_options_create(FC3D_NSGS);\n");
  printf("       int info = fc3d_driver(problem, reaction, velocity, options);\n");
  printf("\n");
}

/* ===========================================================================
 * Convenience Function: Get Solver Info
 * =========================================================================== */

void fc3d_print_solver_info(solver_id_t solver_id) {
  const SolverEntry* solver = solver_registry_lookup(solver_id);
  
  if (!solver) {
    printf("Solver ID %d not found in registry.\n", solver_id);
    printf("Use fc3d_list_available_solvers() to see available solvers.\n");
    return;
  }
  
  printf("\n");
  printf("Solver Information:\n");
  printf("===================\n");
  printf("  Name:        %s\n", solver->name);
  printf("  ID:          %d\n", solver->id);
  printf("  Description: %s\n", solver->description);
  printf("  Type:        %s\n", solver->is_local_solver ? "Local (one-contact)" : "Main solver");
  printf("\n");
  printf("Default Parameters:\n");
  printf("  Max iterations: %d\n", solver->default_max_iter);
  printf("  Tolerance:      %.2e\n", solver->default_tol);
  printf("\n");
  printf("Functions:\n");
  printf("  Init:   %s\n", solver->init ? "Yes" : "No");
  printf("  Solve:  Yes\n");
  printf("  Free:   %s\n", solver->free ? "Yes" : "No");
  printf("  Error:  %s\n", solver->error ? "Yes" : "No");
  printf("\n");
}
