/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef COHESIVE_FRICTION_3D_COMPUTE_ERROR_H
#define COHESIVE_FRICTION_3D_COMPUTE_ERROR_H

/*!\file cohesive_friction_3d_compute_error.h
  \brief Error computation functions for cohesive friction-contact 3D problems

  This file provides functions for computing the error/residual in cohesive
  friction-contact problems, similar to fc3d_compute_error.h for standard
  friction-contact problems.
*/

#include "CohesiveFrictionContactProblem.h"
#include "SolverOptions.h"

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Error computation for a cohesive friction-contact 3D problem
 *
 * Computes the residual using the normal map formulation:
 * error = ||r - P_C(r - u)|| where P_C is projection on friction cone
 *
 * \param[in] problem the cohesive friction-contact problem
 * \param[in] reaction the reaction vector r
 * \param[in] velocity the velocity vector u
 * \param[in] tolerance tolerance for convergence check
 * \param[in] norm normalization factor (typically norm of q)
 * \param[out] error the computed error value
 * \return 0 if error <= tolerance (converged), 1 otherwise
 */
int cohesive_friction_3d_compute_error(CohesiveFrictionContactProblem* problem,
                                       double* reaction, double* velocity, double tolerance,
                                       SolverOptions* options, double norm, double* error);

/**
 * Unitary error computation for a single cohesive friction-contact
 *
 * Computes the error contribution from a single contact using
 * the normal map residual formulation.
 *
 * \param[in] r the local reaction (size 3)
 * \param[in] u the local velocity (size 3)
 * \param[in] mu the friction coefficient
 * \param[in,out] error accumulator for error sum of squares
 * \param[out] worktmp work vector (size 3)
 */
void cohesive_friction_3d_unitary_compute_and_add_error(double r[3], double u[3], double mu,
                                                        double* error, double worktmp[3]);

/**
 * Compute dual cone error for a single cohesive friction-contact
 *
 * Uses the dual cone formulation for error computation.
 *
 * \param[in] r the local reaction (size 3)
 * \param[in] u the local velocity (size 3)
 * \param[in] mu the friction coefficient
 * \param[in,out] error accumulator for error sum of squares
 * \param[out] worktmp work vector (size 3)
 */
void cohesive_friction_3d_unitary_compute_dual_and_add_error(double r[3], double u[3],
                                                             double mu, double* error,
                                                             double worktmp[3]);

#if defined(__cplusplus)
}
#endif

#endif /* COHESIVE_FRICTION_3D_COMPUTE_ERROR_H */
