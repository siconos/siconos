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

#ifndef MIXEDPROBLEM_NSGS_H
#define MIXEDPROBLEM_NSGS_H

/*!\file MixedProblem_NSGS.h
  Generic Non-Smooth Gauss-Seidel (NSGS) solver for MixedProblems.
  Handles heterogeneous combinations of friction, rolling friction,
  and plasticity models.
*/

#include "MixedProblem.h"
#include "MixedProblems_options.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/**
 * Solve a MixedProblem using Non-Smooth Gauss-Seidel (NSGS).
 *
 * This is a generic solver that handles heterogeneous model types
 * (friction 2D/3D, rolling friction 2D/3D, Von Mises, Drucker-Prager)
 * in a single Gauss-Seidel iteration.
 *
 * \param[in] problem the MixedProblem to solve
 * \param[in,out] z the solution vector (size: problem->total_dimension)
 * \param[in,out] w the auxiliary vector (size: problem->total_dimension)
 * \param[out] info return code (0 if successful)
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS(MixedProblem* problem, double* z, double* w, int* info,
                       SolverOptions* options);

/**
 * Initialize the MixedProblem NSGS solver.
 *
 * \param[in] problem the MixedProblem
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS_init(MixedProblem* problem, SolverOptions* options);

/**
 * Free memory allocated for the MixedProblem NSGS solver.
 *
 * \param[in] problem the MixedProblem
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS_free(MixedProblem* problem, SolverOptions* options);

/**
 * Set default solver options for MixedProblem NSGS.
 *
 * \param[in] problem the MixedProblem
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS_set_default_solver_options(MixedProblem* problem,
                                                  SolverOptions* options);

/**
 * Compute error/residual for a MixedProblem.
 *
 * \param[in] problem the MixedProblem
 * \param[in] z the solution vector
 * \param[in] w the auxiliary vector
 * \param[in] tolerance error tolerance
 * \param[in,out] options solver options
 * \param[out] error the computed error
 * \return 0 if successful
 */
int mixedProblem_compute_error(MixedProblem* problem, double* z, double* w, double tolerance,
                               SolverOptions* options, double* error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* MIXEDPROBLEM_NSGS_H */
