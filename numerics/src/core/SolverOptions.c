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

#include "SolverOptions.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for INFINITY
#include <stdio.h>   // for NULL, size_t, printf
#include <stdlib.h>  // for free, calloc, malloc
#include <string.h>  // for strcmp

/* Solver registration system (NEW) */
#include "numerics_verbose.h"
#include "solver_registry.h"

/** Create a struct SolverOptions and initialize its content.

    Allocate memory for internal parameters array, ensure that pointers
    are properly set/nullified and tolerance, max iter and common default values.
    All other specific parameters must be set explicitely for each solver by
    its internal set default function.

    \param options struct to fill
    \param solver_id id number of the solver
    \param iter_max maximum number of iterations allowed for this solver
    \param tol tolerance for the solution.
*/
static SolverOptions* solver_options_initialize(int solver_id, int iter_max, double tol) {
  // Warning : the C structure members after the first instanciation might be anything.
  // (C structure members can not have default values).
  // Maybe we should write a function to create so, like so * = so_new ?
  SolverOptions* options = (SolverOptions*)calloc(1, sizeof(SolverOptions));

  options->solverId = solver_id;
  options->iSize = OPTIONS_PARAM_SIZE;
  options->dSize = OPTIONS_PARAM_SIZE;
  options->iparam = (int*)calloc(options->iSize, sizeof(int));
  options->dparam = (double*)calloc(options->dSize, sizeof(double));

  // The content of iparam and dparam for indices in SICONOS_IPARAM and SICONOS_DPARAM enums
  // enum must be set in this function with a default value.

  options->iparam[SICONOS_IPARAM_MAX_ITER] = iter_max;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->iparam[SICONOS_IPARAM_PREALLOC] = 0;

  options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0;      // not common to all formulations ?
  options->iparam[SICONOS_IPARAM_ERROR_EVALUATION] = 0;  // not common to all formulations ?
  options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] =
      0;  // not common to all formulations ?

  options->dparam[SICONOS_DPARAM_TOL] = tol;
  options->dparam[SICONOS_DPARAM_RESIDU] = INFINITY;

  options->filterOn = true;

  // Everything regarding dwork, iwork, solverData and solverParameters is problem/formulation
  // specific and must be handled later, locally, by the driver when the problem is known.
  //
  options->dWorkSize = 0;
  options->dWork = NULL;
  options->dWorkSize = 0;
  options->iWork = NULL;
  options->callback = NULL;
  options->numberOfInternalSolvers = 0;
  // The number of internal solvers and the required allocation must be done in
  // ..._set_default function, not here
  // options->numberOfInternalSolvers = number_of_internal_solvers;
  // if (options->numberOfInternalSolvers > 0) {
  //   options->internalSolvers =
  //       calloc(options->numberOfInternalSolvers, sizeof(SolverOptions*));
  // }
  options->solverData = NULL;
  options->solverParameters = NULL;

  options->isSet = true;
  return options;
}

/** */
static void recursive_solver_options_print(SolverOptions* options, int level) {
  if (!options) {
    numerics_printf("%s========== NULL options pointer.\n");
    return;
  }
  char* marge;
  marge = (char*)malloc((size_t)(level + 1) * sizeof(char));
  for (int i = 0; i < level; i++) marge[i] = ' ';
  marge[level] = '\0';

  numerics_printf("%s========== solver parameters: ", marge);
  if (!options->isSet)
    numerics_printf("%sThe solver parameters have not been set. \t options->isSet = %i ",
                    marge, options->isSet);
  else {
    numerics_printf("%sThe solver parameters below have  been set \t options->isSet = %i",
                    marge, options->isSet);
    numerics_printf("%sId of the solver\t\t\t\t options->solverId = %d ", marge,
                    options->solverId);
    numerics_printf("%sName of the solver\t\t\t\t %s ", marge,
                    solver_options_id_to_name(options->solverId));
    if (options->iparam != NULL) {
      numerics_printf("%ssize of the int parameters\t\t\t options->iSize = %zu", marge,
                      options->iSize);
      numerics_printf("%snon zero int parameters in options->iparam:", marge);
      for (size_t i = 0; i < options->iSize; ++i) {
        if (options->iparam[i])
          numerics_printf("%s\t\t\t\t\t\t options->iparam[%i] = %d", marge, i,
                          options->iparam[i]);
      }
    }
    if (options->dparam != NULL) {
      numerics_printf("%sdouble parameters \t\t\t\t options->dparam", marge);
      numerics_printf("%ssize of the double parameters\t\t\t options->dSize = %zu", marge,
                      options->dSize);
      numerics_printf("%snon zero double parameters in options->dparam:", marge);
      for (size_t i = 0; i < options->dSize; ++i) {
        if (options->dparam[i] > 0.)
          numerics_printf("%s\t\t\t\t\t\t options->dparam[%i] = %.6le", marge, i,
                          options->dparam[i]);
      }
    }
  }
  if (options->iWork == NULL) {
    numerics_printf("%sinteger work array have not been allocated. \t options->iWork = NULL ",
                    marge);
  } else {
    numerics_printf("%sinteger work array have been allocated. \t options->iWork = %p ", marge,
                    options->iWork);
    numerics_printf("%sinteger work array size \t\t\t options->iSize = %zu ", marge,
                    options->iSize);
  }
  if (options->dWork == NULL) {
    numerics_printf("%sdouble work array have not been allocated. \t options->dWork = NULL ",
                    marge);
  } else {
    numerics_printf("%sdouble work array have been allocated. \t options->dWork = %p ", marge,
                    options->dWork);
    numerics_printf("%sdouble work array size \t\t\t options->dSize = %zu ", marge,
                    options->dSize);
  }

  numerics_printf("%sSee %s documentation for parameters definition)", marge,
                  solver_options_id_to_name(options->solverId));

  numerics_printf(
      "%snumber of internal (or local) solvers \t\t options->numberOfInternalSolvers = %i",
      marge, options->numberOfInternalSolvers);
  for (size_t i = 0; i < options->numberOfInternalSolvers; i++) {
    recursive_solver_options_print(options->internalSolvers[i], level + 1);
  }
  free(marge);
}

void solver_options_print(SolverOptions* options) {
  recursive_solver_options_print(options, 0);
}

void solver_options_delete(SolverOptions* op) {
  if (op) {
    // Clear solverParameters and solverData, before anything.
    // Remark : these are specific data. And so, alloc/release
    // memory operations should be handled inside each
    // solver, at the end of the driver call.
    if (op->solverParameters) free(op->solverParameters);
    op->solverParameters = NULL;

    if (op->solverData) free(op->solverData);
    op->solverData = NULL;

    // Clear callback
    if (op->callback) free(op->callback);
    op->callback = NULL;

    // Clear internal solver(s)
    if (op->internalSolvers) {
      for (size_t i = 0; i < op->numberOfInternalSolvers; i++) {
        solver_options_delete(op->internalSolvers[i]);
        op->internalSolvers[i] = NULL;
      }
      free(op->internalSolvers);
      op->numberOfInternalSolvers = 0;
    }
    op->internalSolvers = NULL;

    // working arrays
    // Solver-specific data : see remark on top of this function,
    // regarding solverData.
    if (op->iWork) free(op->iWork);
    op->iWork = NULL;
    if (op->dWork) free(op->dWork);
    op->dWork = NULL;

    // int and double parameters arrays
    if (op->iparam) free(op->iparam);
    op->iparam = NULL;
    if (op->dparam) free(op->dparam);
    op->dparam = NULL;

    op->isSet = false;
    free(op);
  }
}

SolverOptions* solver_options_copy(SolverOptions* source) {
  // Create a new solver options, with default setup
  SolverOptions* options = solver_options_create(source->solverId);

  for (size_t i = 0; i < OPTIONS_PARAM_SIZE; ++i) {
    options->iparam[i] = source->iparam[i];
    options->dparam[i] = source->dparam[i];
  }

  if (source->dWork) {
    options->dWork = calloc(source->dWorkSize, sizeof(double));
    options->dWorkSize = source->dWorkSize;
    for (size_t i = 0; i < options->dWorkSize; ++i) options->dWork[i] = source->dWork[i];
  }
  if (source->iWork) {
    options->iWork = calloc(source->iWorkSize, sizeof(int));
    options->iWorkSize = source->iWorkSize;
    for (size_t i = 0; i < options->iWorkSize; ++i) options->iWork[i] = source->iWork[i];
  }
  assert(options->numberOfInternalSolvers == source->numberOfInternalSolvers);
  // this assert should be ensured by solver_options_create and initialize.

  for (size_t i = 0; i < options->numberOfInternalSolvers; ++i)
    options->internalSolvers[i] = solver_options_copy(source->internalSolvers[i]);

  // Warning pointer links!
  if (source->callback)
    options->callback =
        source->callback;  // Note FP: is it really safe to create pointer link here?

  if (source->solverData) options->solverData = source->solverData;

  if (source->solverParameters) options->solverParameters = source->solverParameters;

  return options;
}

void solver_options_update_internal(SolverOptions* parent, size_t internal_solver_number,
                                    int solver_id) {
  // Avoid access to a non existing internal solver.
  assert(parent->numberOfInternalSolvers > internal_solver_number);

  // Destroy current internal solver and create a new one with the new id.
  solver_options_delete(parent->internalSolvers[internal_solver_number]);
  parent->internalSolvers[internal_solver_number] = solver_options_create(solver_id);
}

const char* solver_options_id_to_name(int Id) {
  const SolverEntry* solver = solver_registry_lookup(Id);
  if (solver) {
    return solver->name;
  }
  return "NONAME";
}

int solver_options_name_to_id(const char* pName) {
  const SolverEntry* solver = solver_registry_lookup_by_name(pName);
  if (solver) {
    return solver->id;
  }
  return 0;
}

SolverOptions* solver_options_get_internal_solver(SolverOptions* options, size_t n) {
  if (n + 1 > options->numberOfInternalSolvers) {
    printf(
        "solver_options_get_internal_solver : the index must be between 0 and  "
        "options->numberOfInternalSolvers -1 ");
    return NULL;
  } else
    return options->internalSolvers[n];
}

void solver_options_set_internal_solver(SolverOptions* options, size_t n, SolverOptions* NSO) {
  if (n + 1 > options->numberOfInternalSolvers) {
    printf(
        "solver_options_set_internal_solver : the index must be between 0 and  "
        "options->numberOfInternalSolvers -1 ");
  } else {
    options->internalSolvers[n] = NSO;
  }
}

/* ===========================================================================
 * Registry-based Solver Options (NEW)
 * =========================================================================== */

SolverOptions* solver_options_create(int solver_id) {
  /* Look up solver in registry */
  const SolverEntry* solver = solver_registry_lookup(solver_id);
  if (!solver) {
    fprintf(stderr, "solver_options_create: solver '%i' not found in registry : %s\n",
            solver_id, solver_options_id_to_name(solver_id));
    return NULL;
  }

  /* Create options with registered defaults */
  SolverOptions* options =
      solver_options_initialize(solver_id, solver->default_max_iter, solver->default_tol);

  /* Call set_default if provided - this sets solver-specific options
   * and allocates internal solvers without requiring a problem pointer */
  if (solver->set_default) {
    solver->set_default(options);
  }
  /* Note: We don't call solver->init here because many init functions
   * require a valid problem pointer (not NULL). The init is typically
   * called later in the driver before solving.
   */

  return options;
}

SolverOptions* solver_options_create_by_name(const char* solver_name) {
  if (!solver_name) return NULL;

  /* Look up solver by name */
  const SolverEntry* solver = solver_registry_lookup_by_name(solver_name);
  if (!solver) {
    fprintf(stderr, "solver_options_create_by_name: solver '%s' not found in registry\n",
            solver_name);
    return NULL;
  }

  return solver_options_create(solver->id);
}

SolverOptions* solver_options_create_and_init(int solver_id, void* problem) {
  SolverOptions* options = solver_options_create(solver_id);
  if (!options) return NULL;

  /* Call solver init with problem if provided */
  const SolverEntry* solver = solver_registry_lookup(solver_id);
  if (solver && solver->init && problem) {
    int init_status = solver->init(problem, options);
    if (init_status != 0) {
      fprintf(stderr, "solver_options_create_and_init: solver init failed with status %d\n",
              init_status);
      solver_options_delete(options);
      return NULL;
    }
  }

  return options;
}

void solver_options_reset_to_defaults(SolverOptions* options) {
  if (!options) return;

  /* Look up solver in registry */
  const SolverEntry* solver = solver_registry_lookup(options->solverId);
  if (!solver) return;

  /* Reset to registered defaults */
  options->iparam[SICONOS_IPARAM_MAX_ITER] = solver->default_max_iter;
  options->dparam[SICONOS_DPARAM_TOL] = solver->default_tol;
}
