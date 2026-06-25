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
#ifndef COHESIVE_FRICTION_3D_NSGS_H
#define COHESIVE_FRICTION_3D_NSGS_H

/*!\file cohesive_friction_3d_nsgs.h
  Non-Smooth Gauss-Seidel (NSGS) solver for 3D cohesive friction-contact problems.

  The NSGS solver iterates over all contacts, solving a local cohesive
  friction-contact problem at each contact while keeping other reactions fixed.
*/

#include "CohesiveFrictionContactProblem.h"
#include "SolverOptions.h"

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * NSGS solver for 3D cohesive friction-contact problems
 *
 * This solver implements the Non-Smooth Gauss-Seidel method with support
 * for cohesive forces. At each iteration, it solves a local friction-contact
 * problem with cohesion for each contact.
 *
 * \param[in,out] problem the cohesive friction-contact problem
 * \param[out] reaction the reaction vector (solution)
 * \param[out] velocity the velocity vector
 * \param[in,out] options solver options
 * \return 0 if converged, error code otherwise
 *
 * \note The cohesion vectors c_n and c_t in the problem are used to modify 
 *       the effective q vector through matrices V, X, U.
 */
int cohesive_friction_3d_nsgs(CohesiveFrictionContactProblem *problem,
                              double *reaction,
                              double *velocity,
                              SolverOptions *options);


#if defined(__cplusplus)
}
#endif

#endif  // COHESIVE_FRICTION_3D_NSGS_H
