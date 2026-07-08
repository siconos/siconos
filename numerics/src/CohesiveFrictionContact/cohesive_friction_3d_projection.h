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
#ifndef COHESIVE_FRICTION_3D_PROJECTION_H
#define COHESIVE_FRICTION_3D_PROJECTION_H

/*!\file cohesive_friction_3d_projection.h
  Projection-based solvers for 3D cohesive friction-contact problems.

  These solvers use projection onto the friction cone to solve the
  cohesive friction-contact problem. They can be used as local solvers
  within NSGS or as standalone solvers.
*/

#include "CohesiveFrictionContactProblem.h"
#include "SolverOptions.h"

#if defined(__cplusplus)
extern "C" {
#endif

void cohesive_friction_3d_projection_initialize(CohesiveFrictionContactProblem* main_problem,
                                                SolverOptions* localsolver_option);

void cohesive_friction_3d_projection_free(CohesiveFrictionContactProblem* main_problem,
                                          SolverOptions*);

/**
 * Project a reaction onto the friction cone
 *
 * \param[in,out] reaction the reaction vector (input/output)
 * \param[in] mu the friction coefficient
 * \param[in] dim the dimension (2 or 3)
 *
 */
int cohesive_friction_3d_projection_solve(CohesiveFrictionContactProblem* localproblem,
                                          double* reaction, SolverOptions* options);
/**
 * Projection on cone solver for 3D cohesive friction-contact
 *
 * This solver uses a fixed-point iteration with projection onto the
 * friction cone. It includes the cohesive force contribution.
 *
 * \param[in,out] problem the cohesive friction-contact problem
 * \param[out] reaction the reaction vector (solution)
 * \param[out] velocity the velocity vector
 * \param[in,out] options solver options
 * \return 0 if converged, error code otherwise
 */

/**
 * Projection on cone with local iteration
 *
 * This solver uses a fixed-point iteration with projection, but performs
 * additional local iterations for each contact to improve convergence.
 *
 * \param[in,out] problem the cohesive friction-contact problem
 * \param[out] reaction the reaction vector (solution)
 * \param[out] velocity the velocity vector
 * \param[in,out] options solver options
 * \return 0 if converged, error code otherwise
 */

#if defined(__cplusplus)
}
#endif

#endif  // COHESIVE_FRICTION_3D_PROJECTION_H
