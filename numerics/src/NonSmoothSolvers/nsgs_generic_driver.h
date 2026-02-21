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

/*!\file nsgs_generic_driver.h
 * \brief Driver for using generic NSGS with different problem types
 *
 * This header provides a unified interface for calling the generic NSGS
 * solver from different problem-specific drivers.
 */

#ifndef NSGS_GENERIC_DRIVER_H
#define NSGS_GENERIC_DRIVER_H

#include "nsgs_generic.h"

/** Generic NSGS solver function signature */
typedef void (*NSGSSolver)(void* problem, double* reaction, double* velocity, int* info,
                           SolverOptions* options);

/** Structure to hold problem-specific callbacks */
typedef struct {
  const char* problem_name;          /* Problem type identifier */
  int dimension;                     /* Dimension per contact */
  NSGSUpdateLocalProblem update;     /* Update local problem */
  NSGSSolveLocal solve;              /* Solve local problem */
  NSGSComputeError compute_error;    /* Compute error */
  NSGSAllocLocalProblem alloc_local; /* Allocate local problem */
  NSGSAllocShuffled alloc_shuffled;  /* Allocate shuffled contacts (optional) */
  NSGSAllocFreezing alloc_freezing;  /* Allocate freezing contacts (optional) */
  NSGSSetLocalTol set_tol;           /* Set local tolerance (optional) */
  NSGSIncrError incremental_error;   /* Incremental error function (optional) */
  NSGSAcceptLocal accept_local;      /* Accept local reaction (optional) */
  NSGSStatsCallback stats;           /* Statistics callback (optional) */
} NSGSProblemCallbacks;

/** Generic NSGS driver
 *
 * This is the main entry point for using NSGS with any problem type.
 * It configures the toolkit from callbacks and calls the generic solver.
 *
 * \param[in] problem Problem-specific structure
 * \param[in,out] reaction Global reaction vector
 * \param[in,out] velocity Global velocity vector
 * \param[out] info Solver status
 * \param[in] options Solver options
 * \param[in] callbacks Problem-specific callbacks
 * \param[in] nc Number of contacts
 * \param[in] q Right-hand side vector
 * \param[in] M Matrix
 */
static inline void nsgs_driver(void* problem, double* reaction, double* velocity, int* info,
                               SolverOptions* options, NSGSProblemCallbacks* callbacks,
                               unsigned int nc, double* q, void* M) {
  /* Setup problem data */
  NSGSProblemData problem_data = {.nc = nc,
                                  .q = q,
                                  .M = M,
                                  .mu = NULL,
                                  .mu_r = NULL,
                                  .storage_type = 0,
                                  .dimension = callbacks->dimension};

  /* Setup toolkit from callbacks */
  NSGSLocalToolkit toolkit = {.update_local_problem = callbacks->update,
                              .solve_local = callbacks->solve,
                              .compute_error = callbacks->compute_error,
                              .incremental_error = callbacks->incremental_error,
                              .accept_local = callbacks->accept_local,
                              .check_convergence = NULL,
                              .alloc_freezing = callbacks->alloc_freezing,
                              .alloc_shuffled = callbacks->alloc_shuffled,
                              .alloc_local = callbacks->alloc_local,
                              .set_local_tol = callbacks->set_tol,
                              .stats_callback = callbacks->stats,
                              .dimension = callbacks->dimension,
                              .omega = 1.0,
                              .use_freezing = 0,
                              .use_shuffling = 0,
                              .use_incremental_error = (callbacks->incremental_error != NULL),
                              .filter_local_sol = 1,
                              .error_eval_type = 0};

  /* Extract options from SolverOptions if needed */
  if (options) {
    int* iparam = options->iparam;
    double* dparam = options->dparam;

    /* Check for common option names across problem types */
    if (iparam) {
      /* These indices should be standardized across problems */
      toolkit.use_freezing = iparam[7];  /* SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT */
      toolkit.use_shuffling = iparam[5]; /* SICONOS_FRICTION_3D_NSGS_SHUFFLE */
      toolkit.filter_local_sol =
          iparam[8];                       /* SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION */
      toolkit.error_eval_type = iparam[9]; /* SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION */
    }

    if (dparam) {
      toolkit.omega = dparam[10]; /* SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE */
    }
  }

  /* Call generic solver */
  nsgs_solve(problem, reaction, velocity, info, options, &toolkit, &problem_data);
}

#endif /* NSGS_GENERIC_DRIVER_H */
