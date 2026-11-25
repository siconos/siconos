// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
int cohesive_friction_3d_nsgs(CohesiveFrictionContactProblem* problem, double* reaction,
                              double* velocity, SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif  // COHESIVE_FRICTION_3D_NSGS_H
