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

/*!\file rolling_fc3d_driver.c
 * \brief RFC3D driver using the solver registration system
 */

#include "RollingFrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "rolling_fc_Solvers.h"
#include "rfc3d_short_names.h"
#include "numerics_verbose.h"

/* Registration-based headers */
#include "solver_registry.h"
#include "numerics_errors.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* ===========================================================================
 * Trivial Case Check
 * =========================================================================== */

int rolling_fc3d_checkTrivialCase(RollingFrictionContactProblem* problem, double* velocity,
                                  double* reaction, SolverOptions* options) {
  (void)options;
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  int n = 5 * nc;
  
  for (int i = 0; i < nc; i++) {
    if (q[5 * i] < -DBL_EPSILON) return -1;
  }
  for (int i = 0; i < n; ++i) {
    velocity[i] = q[i];
    reaction[i] = 0.;
  }
  
  numerics_printf("rolling_fc3d: trivial solution (take-off), reaction = 0, velocity = q.");
  return 0;
}

/* ===========================================================================
 * Registration-Based Driver
 * =========================================================================== */

int rolling_fc3d_driver(RollingFrictionContactProblem* problem, double* reaction,
                        double* velocity, SolverOptions* options) {
  /* Input validation using standardized macros */
  CHECK_NULL(problem);
  CHECK_NULL(reaction);
  CHECK_NULL(velocity);
  CHECK_OPTIONS(options);
  CHECK_NULL(problem->M);
  CHECK_NULL(problem->q);
  CHECK_NULL(problem->mu);
  CHECK_ARG(problem->numberOfContacts > 0, "Number of contacts must be positive");

  /* Check dimension */
  CHECK_DIMENSION(problem->dimension, 5);

  /* Initialize output */
  SET_SOLVER_ITER_DONE(options, 0);
  SET_SOLVER_RESIDUAL(options, 0.0);

  /* Check for trivial case */
  int trivial_status = rolling_fc3d_checkTrivialCase(problem, velocity, reaction, options);
  if (trivial_status == 0) {
    return NUMERICS_OK;
  }

  /* Lookup solver in registry */
  const SolverEntry* solver = solver_registry_lookup(options->solverId);

  if (!solver) {
    numerics_printf("rolling_fc3d_driver: solver ID %d not found in registry", options->solverId);
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  numerics_printf_verbose(1, "rolling_fc3d_driver: using solver '%s' (%s)",
                          solver->name, solver->description);

  /* Validate solver is appropriate for this problem type */
  if (solver->is_local_solver) {
    numerics_printf("rolling_fc3d_driver: solver '%s' is a local solver, cannot be used as main solver",
                    solver->name);
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  /* Check solve function exists */
  if (!solver->solve) {
    numerics_printf("rolling_fc3d_driver: solver '%s' has no solve function", solver->name);
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  /* Initialize solver if init function provided */
  if (solver->init) {
    int init_status = solver->init(problem, options);
    if (init_status != NUMERICS_OK) {
      numerics_printf("rolling_fc3d_driver: solver initialization failed with error %d", init_status);
      return init_status;
    }
  }

  /* Call the solver */
  numerics_printf_verbose(1, "rolling_fc3d_driver: calling solver...");
  int solve_status = solver->solve(problem, reaction, velocity, options);

  /* Cleanup if needed */
  if (solver->free) {
    solver->free(problem, options);
  }

  /* Log result */
  if (solve_status == NUMERICS_OK) {
    numerics_printf_verbose(1, "rolling_fc3d_driver: solver converged successfully");
  } else {
    numerics_printf_verbose(1, "rolling_fc3d_driver: solver returned status %d", solve_status);
  }

  return solve_status;
}
