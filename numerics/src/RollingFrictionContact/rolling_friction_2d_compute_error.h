/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef rolling_friction_3d_compute_error_H
#define rolling_friction_3d_compute_error_H

/*!\file rolling_friction_3d_compute_error.h
  \brief functions related to error computation for friction-contact problems

*/

#include "NumericsFwd.h"  // for RollingFrictionContactProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/** Error computation (using the normal map residual) for rolling friction-contact 2D problem
    \param problem the structure which defines the friction-contact problem
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param norm norm of a vector (problem->q) for relative error
    \param[in,out] error value
    \return 0 if ok
 */
int rolling_friction_2d_compute_error(RollingFrictionContactProblem* problem, double* z,
                                      double* w, double tolerance, SolverOptions* options,
                                      double norm, double* error);

/** Error computation (using the normal map residual) for one rolling friction-contact 2D
   problem \param r the reaction force \param u the local velocity \param mu coeficient of
   friction \param worktmp work vector \param[in,out] error value
 */
void rolling_friction_2d_unitary_compute_and_add_error(double r[3], double u[3], double mu,
                                                       double mur, double* error,
                                                       double* worktmp);

#if defined(__cplusplus)
}
#endif

#endif
