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

/*!\file rolling_friction_2d_driver.c
 * \brief RFC2D driver using the solver registration system
 */

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "naming_conventions.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "rolling_fc_Solvers.h"
#include "solver_registry.h"

/* ===========================================================================
 * Trivial Case Check
 * =========================================================================== */

int rolling_friction_2d_checkTrivialCase(RollingFrictionContactProblem* problem,
                                         double* velocity, double* reaction,
                                         SolverOptions* options) {
  (void)options;
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  int n = 3 * nc;

  for (int i = 0; i < nc; i++) {
    if (q[3 * i] < -DBL_EPSILON) return -1;
  }
  for (int i = 0; i < n; ++i) {
    velocity[i] = q[i];
    reaction[i] = 0.;
  }

  numerics_printf("rolling_fc2d: trivial solution (take-off), reaction = 0, velocity = q.");
  return 0;
}

/* ===========================================================================
 * Registration-Based Driver
 * =========================================================================== */

int rolling_friction_2d_driver(RollingFrictionContactProblem* problem, double* reaction,
                               double* velocity, SolverOptions* options) {
  /* Input validation */
  if (!problem || !reaction || !velocity || !options) {
    numerics_error_log("rolling_friction_2d_driver", "null input argument");
    return NUMERICS_ERR_NULL_POINTER;
  }

  /* Check dimension */
  CHECK_DIMENSION(problem->dimension, 3);

  /* Initialize output */
  SET_SOLVER_ITER_DONE(options, 0);
  SET_SOLVER_RESIDUAL(options, 0.0);

  /* Check for trivial case */
  int trivial_status =
      rolling_friction_2d_checkTrivialCase(problem, velocity, reaction, options);
  if (trivial_status == 0) {
    return NUMERICS_OK;
  }

  /* Lookup solver in registry */
  const SolverEntry* solver = solver_registry_lookup(options->solverId);

  CHECK_COND(solver != NULL, NUMERICS_ERR_INVALID_SOLVER, "Solver not found");

  numerics_printf_verbose(1, "rolling_friction_2d_driver: using solver '%s' (%s)",
                          solver->name, solver->description);

  /* Validate solver is appropriate for this problem type */
  if (solver->is_local_solver) {
    numerics_printf(
        "rolling_friction_2d_driver: solver '%s' is a local solver, cannot be used as main "
        "solver",
        solver->name);
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  /* Check solve function exists */
  CHECK_COND(solver->solve != NULL, NUMERICS_ERR_INVALID_SOLVER,
             "Solver has no solve function");

  /* Initialize solver if init function provided */
  if (solver->init) {
    int init_status = solver->init(problem, options);
    if (init_status != NUMERICS_OK) {
      numerics_printf("rolling_friction_2d_driver: solver initialization failed with error %d",
                      init_status);
      return init_status;
    }
  }

  /* Call the solver */
  numerics_printf_verbose(1, "rolling_friction_2d_driver: calling solver...");
  int solve_status = solver->solve(problem, reaction, velocity, options);

  /* Cleanup if needed */
  if (solver->free) {
    solver->free(problem, options);
  }

  /* Log result */
  if (solve_status == NUMERICS_OK) {
    numerics_printf_verbose(1, "rolling_friction_2d_driver: solver converged successfully");
  } else {
    numerics_printf_verbose(1, "rolling_friction_2d_driver: solver returned status %d",
                            solve_status);
  }

  return solve_status;
}
