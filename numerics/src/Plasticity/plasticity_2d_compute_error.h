// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef plasticity_2d_compute_error_H
#define plasticity_2d_compute_error_H

/*!\file plasticity_2d_compute_error.h
  \brief functions related to error computation for Mohr Coulomb 2D plasticity

*/

#include "NumericsFwd.h"  // for PlasticityProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/**
    Error computation (using the normal map residual) for Mohr Coulomb 2D plasticity

    \param problem the structure which defines the  Mohr Coulomb 2D plasticity
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param norm norm of a vector (problem->q) for relative error
    \param[in,out] error value
    \return 0 if ok
 */
int plasticity_2d_compute_error(PlasticityProblem* problem, double* z, double* w,
                                double tolerance, SolverOptions* options, double norm,
                                double* error);

/**
    Error computation (using the normal map residual) for one Mohr Coulomb 2D plasticity

    \param r the reaction force
    \param u the local velocity
    \param eta coefficient of the cone
    \param theta coefficient of the dilatency
    \param worktmp work vector
    \param[in,out] error value
 */
void plasticity_2d_unitary_compute_and_add_error(double* restrict r, double* restrict u,
                                                 double eta, double theta, double* error,
                                                 double* worktmp);

/** Check whether sigma = 0 solves the one-cone Drucker-Prager problem with w = q:
 *  q + theta ||q_T|| e_N must be in the dual cone K_eta^* = {eta ||v_T|| <= v_N}, i.e.
 *  eta ||q_T|| <= q_N + theta ||q_T||. (q_N >= 0 is enough only when theta >= eta.)
 *  \param[in] q the local vector (w for sigma = 0)
 *  \param[in] eta cone coefficient
 *  \param[in] theta dilatancy coefficient
 *  \return 1 if sigma = 0 is a solution, 0 otherwise
 */
int plasticity_2d_DruckerPrager_zero_is_solution(const double q[3], double eta, double theta);

#if defined(__cplusplus)
}
#endif

#endif
