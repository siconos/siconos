/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef global_rolling_friction_3d_compute_error_H
#define global_rolling_friction_3d_compute_error_H

/*!\file global_rolling_friction_3d_compute_error.h
  \brief functions related to error computation for rolling friction-contact problems

*/

#include "NumericsFwd.h"  // for GlobalRollingFrictionContactProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/** Error computation (using the normal map residual) for rolling friction-contact 3D problem
    \param problem the structure which defines the rolling friction-contact problem
    \param u vector
    \param r vector
    \param tolerance value for error computation
    \param options
    \param norm of a vector (problem->q) for relative error
    \param[in,out] error value
    \return 0 if ok
 */
int global_rolling_friction_3d_compute_error(GlobalRollingFrictionContactProblem* problem,
                                             double* reaction, double* velocity,
                                             double* globalVelocity, double tolerance,
                                             double* error, int problemIsNotConvex);

/** Error computation (using the normal map residual) for one rolling friction-contact 3D
   problem \param r the reaction force \param u the local velocity \param mu coeficient of
   friction \param mur coeficient of rolling \param worktmp work vector \param[in,out] error
   value
 */
void global_rolling_friction_3d_unitary_compute_and_add_error(double* r, double* u, double mu,
                                                              double mur, double* error,
                                                              double* worktmp,
                                                              int problemIsNotConvex);

#if defined(__cplusplus)
}
#endif

#endif
