/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef VariationalInequality_compute_error_H
#define VariationalInequality_compute_error_H

/*!\file VariationalInequality_computeError.h
  \brief functions related to error computation for friction-contact problems

*/

#include "NumericsFwd.h"  // for VariationalInequality, SolverOptions

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/** Error computation for a VI problem. This function requires dWork to point to
 * at least 2*n double of allocated memory or it malloc this memory
    \param problem the structure which defines the VI problem
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options solver options
    \param[in,out] error value
    \return 0 if ok
 */
int variationalInequality_computeError(VariationalInequality* problem, double* z, double* w,
                                       double tolerance, SolverOptions* options,
                                       double* error);

/** Error computation for a box VI problem, that is\f$ \Pi_box(x-F(x)) - x\f$
    \param problem the structure which defines the VI problem
    \param x vector
    \param F vector
    \param tolerance value for error computation
    \param[in,out] error value
    \return 0 if ok
 */
int variationalInequality_compute_error_box(VariationalInequality* problem, double* restrict x,
                                            double* restrict F, double tolerance,
                                            double* restrict error);

#if defined(__cplusplus)
}
#endif

#endif
