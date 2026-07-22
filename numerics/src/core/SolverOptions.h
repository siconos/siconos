/* Siconos is a program dedicated to modeling, simulation and
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
#ifndef SolverOptions_H
#define SolverOptions_H

/*!\file SolverOptions.h
  Structure used to send options (name, parameters and so on) to a specific
  solver-driver (mainly from Kernel to Numerics).
*/
#include <stdbool.h>  // for boolean type
#include <stdio.h>    // for size_t

#include "NumericsFwd.h"  // for SolverOptions

/* Forward declaration for solver registry (avoid circular include) */
typedef int solver_id_t;

/**
    Structure used to store user callbacks inside solvers
*/
typedef struct {
  void *env; /**< general user environment */
  void (*collectStatsIteration)(
      void *env, size_t size, double *reaction, double *velocity, double error,
      void *extra_data); /**<pointer on a function
                          * Its signature is: user env, problem size, reaction,
                          * velocity, error at end of solver iteration (when
                          * this makes sense) and an extra data structure */
} Callback;

// length of iparam/dparam arrays in solver options.
#define OPTIONS_PARAM_SIZE 30

/**
  Structure used to send options (name, parameters and so on) to a specific
  solver (mainly from Kernel to Numerics).

  Creation, update and destruction:

  - solver_options_create()

  - solver_options_update_internal()

  - solver_options_delete()


  Details in users'guide.

*/
struct SolverOptions {
  int solverId;     /**< id number of the solver. */
  bool isSet;       /**< true(1) if the structure is ready to be used by a numerics
                       driver. */
  size_t iSize;     /**< iSize size of vector iparam */
  int *iparam;      /**< list of solver parameters (integer type); Check solvers doc
                       for details. */
  size_t dSize;     /**< size of vector dparam */
  double *dparam;   /**< list of solver parameters (double type); Check solvers
                       doc for details. */
  bool filterOn;    /**< if true (1), check solution validity after the driver
                       call. Default = 1.    For example if filterOn = 1 for a LCP,
                       lcp_compute_error()    will be called at the end of the
                       process). */
  size_t dWorkSize; /**< size of double type internal work array.*/
  double *dWork;    /**< internal (double type) work array.*/
  size_t iWorkSize; /**< size of integer type internal work array.*/
  int *iWork;       /**< internal (integer type) work array.*/
  size_t numberOfInternalSolvers;         /**< the number of internal or local
                                             'sub-solvers' used by the solver         (size of
                                             internalSolvers) .*/
  struct SolverOptions **internalSolvers; /**< list of internal solver options*/
  Callback *callback;                     /**< pointer to user-defined callback*/
  void *solverParameters;                 /**< additional parameters specific to the solver
                                             (GAMS and NewtonMethod only) */
  void *solverData;                       /**< additional data specific to the solver */
};

/** Some value for iparam index */
enum SICONOS_IPARAM {
  SICONOS_IPARAM_MAX_ITER = 0,
  SICONOS_IPARAM_ITER_DONE = 1,
  SICONOS_IPARAM_PREALLOC = 2,
  SICONOS_IPARAM_NSGS_SHUFFLE = 5,
  SICONOS_IPARAM_ERROR_EVALUATION = 3,  // CHECK IF THERE ARE NO CONFLICT WITH THIS !!
  SICONOS_IPARAM_PATHSEARCH_STACKSIZE = 19
};

/** allowed values for iparam[SICONOS_IPARAM_ERROR_EVALUATION */
enum SICONOS_IPARAM_ERROR_EVALUATION_ENUM {
  /** Complete error computation, including v computation*/
  SICONOS_ERROR_FULL_EVALUATION = 0,
  /** Light error computation with incremental values on r verification of
     absolute error at the end */
  SICONOS_ERROR_LIGHT_EVALUATION = 1,
  /**  only light error computation, do not update v unknown) */
  SICONOS_ERROR_LIGHT_EVALUATION_NO_UPDATE = 2
};

/** Some values for dparam index */
enum SICONOS_DPARAM {
  SICONOS_DPARAM_TOL = 0,
  SICONOS_DPARAM_RESIDU = 1,
  SICONOS_DPARAM_TIME_BEFORE_LOOP = 20, // time before while loop in fc2d_nsgs_graph_permut
  SICONOS_DPARAM_TIME_IN_LOOP = 21, // time inside while loop in fc2d_nsgs_graph_permut
  SICONOS_DPARAM_TIME_AFTER_LOOP = 22 // time after while loop in fc2d_nsgs_graph_permut
};

#if defined(__cplusplus)
extern "C" {
#endif

/**
    screen display of solver parameters

    \param options the structure to be displayed
*/
void solver_options_print(SolverOptions *options);

/**
   Clear and free all pointer members of the structure, then
   release memory

   \param options the structure to be cleared.
*/
void solver_options_delete(SolverOptions *options);

/**
   Copy an existing set of options, to create a new one. Warning : callback,
   solverData and solverParameters of the new structure are pointer links to
   those of the original one!

   \param source an existing solver options structure
   \return a pointer to options set, ready to use by a driver.
*/
SolverOptions *solver_options_copy(SolverOptions *source);

/**
   Change one of the internal solver of a previously defined SolverOptions set.
   Allocate internal memories and set default values for the internal solver.
   Warning : the actual internal solver in position internal_solver_number and
   all its content will be destroyed and replaced by a new one.

   \param parent the top-level SolverOptions which contains the internal solver
   to be updated
   \param internal_solver_number number of the internal solver to
   be update (warning : this is the position in the list of internal solvers,
   not the id!)
   \param solver_id id number of the new internal solver to be
   created/updated
*/
void solver_options_update_internal(SolverOptions *parent, size_t internal_solver_number,
                                    int solver_id);

/** return the id of a solver based on its name
 *
 *  \param pName the name of the solver
 *  \return the id of the solver or 0 if it failed
 */
int solver_options_name_to_id(const char *pName);

/** return the name of a solver given its id
 *
 *  \param Id the id of the solver
 *  \return the name of the solver
 */
const char *solver_options_id_to_name(int Id);

/**
    return the internal solver options set

    \param options parent options
    \param number of the targeted solver
    \return a pointer to the internal solver options set
*/
SolverOptions *solver_options_get_internal_solver(SolverOptions *options, size_t n);

/**
   set internal solver

   \param options parent options
   \param number of the targeted solver
   \param the solver options to be used as internal solver number n
*/
void solver_options_set_internal_solver(SolverOptions *options, size_t n, SolverOptions *NSO);

/* ===========================================================================
 * Registry-based Solver Options (NEW)
 * =========================================================================== */

/**
 * Create solver options using registration system.
 *
 * This uses the solver registry to:
 * 1. Look up solver metadata (name, defaults)
 * 2. Set default parameters from registration
 * 3. Call solver-specific initialization if available
 *
 * \param solver_id The solver ID (e.g., SICONOS_FRICTION_3D_NSGS)
 * \return Pointer to created options, or NULL if solver not found
 */
SolverOptions *solver_options_create(int solver_id);

/**
 * Create solver options by name.
 *
 * Convenience function to create options using solver name instead of ID.
 *
 * \param solver_name The solver name (e.g., "FC3D_NSGS")
 * \return Pointer to created options, or NULL if name not found
 */
SolverOptions *solver_options_create_by_name(const char *solver_name);

/**
 * Create options and apply solver initialization.
 *
 * Creates options AND calls the solver's init function (if provided).
 *
 * \param solver_id The solver ID
 * \param problem Optional problem pointer (for problem-specific init)
 * \return Pointer to created and initialized options, or NULL on failure
 */
SolverOptions *solver_options_create_and_init(int solver_id, void *problem);

/**
 * Reset options to registered defaults.
 *
 * Restores all parameters to their registered default values.
 *
 * \param options The options to reset (modified in place)
 */
void solver_options_reset_to_defaults(SolverOptions *options);

#if defined(__cplusplus)
}
#endif

#endif
