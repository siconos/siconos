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

/*!\file generic_driver_template.h
 * \brief Template for creating registration-based drivers
 *
 * This header provides a macro-based template for creating drivers
 * that use the solver registration system. This eliminates the need
 * to write repetitive driver code for each problem type.
 *
 * Usage example for LCP:
 *   GENERIC_DRIVER_IMPLEMENTATION(lcp, LinearComplementarityProblem, z, w)
 *
 * This creates:
 *   - lcp_driver_registered()
 *   - lcp_solver_options_create()
 *   - lcp_list_available_solvers()
 *   - lcp_print_solver_info()
 */

#ifndef GENERIC_DRIVER_TEMPLATE_H
#define GENERIC_DRIVER_TEMPLATE_H

#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"
#include "numerics_verbose.h"
#include <stdio.h>
#include <stdlib.h>

/* ===========================================================================
 * Macro: GENERIC_DRIVER_IMPLEMENTATION
 * ===========================================================================
 * Creates a complete driver implementation for a problem type.
 *
 * Parameters:
 *   prefix      - Function prefix (e.g., "lcp", "fc3d", "vi")
 *   ProblemType - The problem struct type
 *   sol         - Solution variable name (e.g., "z", "reaction")
 *   aux         - Auxiliary variable name (e.g., "w", "velocity")
 */

#define GENERIC_DRIVER_IMPLEMENTATION(prefix, ProblemType, sol, aux) \
\
int prefix##_driver_registered(ProblemType* problem, double* sol, \
                               double* aux, SolverOptions* options) { \
  NUMERICS_CHECK_NULL(problem); \
  NUMERICS_CHECK_NULL(sol); \
  NUMERICS_CHECK_NULL(aux); \
  NUMERICS_CHECK_NULL(options); \
  \
  SET_SOLVER_ITER_DONE(options, 0); \
  SET_SOLVER_RESIDUAL(options, 0.0); \
  \
  const SolverEntry* solver = solver_registry_lookup(options->solverId); \
  if (!solver) { \
    numerics_printf_verbose(0, #prefix "_driver_registered: solver ID %d not found", \
                            options->solverId); \
    return NUMERICS_ERR_INVALID_SOLVER; \
  } \
  \
  numerics_printf_verbose(1, #prefix "_driver_registered: using solver '%s' (%s)", \
                          solver->name, solver->description); \
  \
  if (solver->is_local_solver) { \
    numerics_printf_verbose(0, #prefix "_driver_registered: '%s' is a local solver", \
                            solver->name); \
    return NUMERICS_ERR_INVALID_SOLVER; \
  } \
  \
  if (solver->init) { \
    int status = solver->init(problem, options); \
    if (status != NUMERICS_OK) return status; \
  } \
  \
  int status = solver->solve(problem, sol, aux, options); \
  \
  if (solver->free) solver->free(problem, options); \
  \
  return status; \
} \
\
SolverOptions* prefix##_solver_options_create(solver_id_t solver_id) { \
  const SolverEntry* solver = solver_registry_lookup(solver_id); \
  if (!solver || solver->is_local_solver) { \
    fprintf(stderr, #prefix "_solver_options_create: invalid solver ID %d\n", solver_id); \
    return NULL; \
  } \
  SolverOptions* options = solver_options_create(solver_id); \
  if (!options) return NULL; \
  if (solver->init) solver->init(NULL, options); \
  return options; \
} \
\
void prefix##_list_available_solvers(void) { \
  printf("\nAvailable " #prefix " solvers:\n"); \
  printf("------------------------\n"); \
  size_t count; \
  const SolverEntry** solvers = solver_registry_get_all(&count); \
  for (size_t i = 0; i < count; i++) { \
    const SolverEntry* s = solvers[i]; \
    printf("  %-20s (ID=%d, max_iter=%d, tol=%.2e)\n", \
           s->name, s->id, s->default_max_iter, s->default_tol); \
  } \
  printf("\n"); \
}

/* ===========================================================================
 * Usage Examples:
 * ===========================================================================
 *
 * For LCP:
 *   GENERIC_DRIVER_IMPLEMENTATION(lcp, LinearComplementarityProblem, z, w)
 *
 * For MLCP:
 *   GENERIC_DRIVER_IMPLEMENTATION(mlcp, MixedLinearComplementarityProblem, z, w)
 *
 * For NCP:
 *   GENERIC_DRIVER_IMPLEMENTATION(ncp, NonlinearComplementarityProblem, z, F)
 *
 * For VI:
 *   GENERIC_DRIVER_IMPLEMENTATION(vi, VariationalInequality, x, F)
 *
 * For QP:
 *   GENERIC_DRIVER_IMPLEMENTATION(qp, ConvexQP, x, F)
 */

#endif /* GENERIC_DRIVER_TEMPLATE_H */
