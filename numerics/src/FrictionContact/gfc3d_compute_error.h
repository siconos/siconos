// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef gfc3d_compute_error_H
#define gfc3d_compute_error_H
#include "NumericsFwd.h"  // for GlobalFrictionContactProblem, SolverOptions

/*!\file gfc3d_compute_error.h
  \brief functions related to error computation for friction-contact problems
*/

#if defined(__cplusplus)
extern "C" {
#endif

/** Error computation for global friction-contact 3D problem
 *
 *  The computation of the error uses as input the reaction (reaction) and the global velocity
 * (globalVelocity) The value of the local velocity (velocity) is recomputed
 *
 *  \param problem the structure which defines the friction-contact problem
 *  \param[in] reaction
 *  \param[in] velocity
 *  \param[out] globalVelocity
 *  \param tolerance value for error computation
 *  \param options pointer to SolverOptions
 *  \param norm_q norm of q or a normalization value
 *  \param norm_b norm of b or a normalization value
 *  \param[in,out] error value
 *  \return 0 if successfull
 */
int gfc3d_compute_error(GlobalFrictionContactProblem* problem, double* reaction,
                        double* velocity, double* globalVelocity, double tolerance,
                        SolverOptions* options, double norm_q, double norm_b, double* error);

int gfc3d_compute_error_convex(GlobalFrictionContactProblem* problem, double* reaction,
                               double* velocity, double* globalVelocity, double tolerance,
                               SolverOptions* options, double norm_q, double norm_b,
                               double* error);

/** Error computation using the norm infinity conic for global friction-contact 3D problem

 *  The computation of the error uses as input the reaction (reaction) and the global velocity
 (globalVelocity)
 *  The value of the local velocity (velocity) is recomputed

 *  \param problem the structure which defines the friction-contact problem
 *  \param[in] reaction
 *  \param[in] velocity
 *  \param[out] globalVelocity
 *  \param tolerance value for error computation
 *  \param options pointer to SolverOptions
 *  \param norm_q norm of q or a normalization value
 *  \param norm_b norm of b or a normalization value
 *  \param[in,out] error value
 *  \return 0 if successfull
 */
int gfc3d_compute_error_norm_infinity_conic(GlobalFrictionContactProblem* problem,
                                            double* reaction, double* velocity,
                                            double* globalVelocity, double tolerance,
                                            SolverOptions* options, double norm_q,
                                            double norm_b, double* error);

/** Function-pointer type for GFC3D error computations.
 *
 *  This allows solvers to plug in different GFC3D error computations.
 */
typedef int (*GFC3D_compute_error_function)(GlobalFrictionContactProblem* problem,
                                            double* reaction, double* velocity,
                                            double* globalVelocity, double tolerance,
                                            SolverOptions* options, double norm_q,
                                            double norm_b, double* error);

#if defined(__cplusplus)
}
#endif

#endif
