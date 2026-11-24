// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef SOCLCP_compute_error_H
#define SOCLCP_compute_error_H

/*!\file soclcp_compute_error.h
  \brief functions related to error computation for SOCLCP

*/
#include "NumericsFwd.h"  // for SecondOrderConeLinearComplementarityProblem

#if defined(__cplusplus)
extern "C" {
#endif

/** Error computation for SOCLCP problem
    \param problem the structure which defines the SOCLCP
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param[in,out] error value
    \return 0 if ok
 */
int soclcp_compute_error(SecondOrderConeLinearComplementarityProblem* problem, double* z,
                         double* w, double tolerance, SolverOptions* options, double* error);

/** Error computation for one SOCLCP problem
    \param z vector
    \param w vector
    \param dim dimension of the cone
    \param mu coeficient of friction
    \param[in,out] error value
    \param worktmp
 */
void soclcp_unitary_compute_and_add_error(double z[3], double w[3], unsigned int dim,
                                          double mu, double* error, double* worktmp);

/** Error computation for SOCLCP problem
    \param problem the structure which defines the SOCLCP
    \param z vector
    \param w vector
    \param options
    \param tolerance value for error computation
    \param[in,out] error value
    \return 0 if ok
 */
int soclcp_compute_error_v(SecondOrderConeLinearComplementarityProblem* problem, double* z,
                           double* w, double tolerance, SolverOptions* options, double* error);

#if defined(__cplusplus)
}
#endif

#endif
