// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef GFC3D_IPM_INTERNAL_H
#define GFC3D_IPM_INTERNAL_H

#include "NumericsMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Returns an array of step-lengths to the boundary reduced by a factor gamma. Uses long
 * double.
 */
double* array_getStepLength(const double* const x, const double* const dx,
                            const size_t vecSize, const size_t varsCount, const double gamma);

/* Returns the primal constraint vector with smoothing variable s for global fricprob:
 * out = velocity - H x globalVelocity - w - phi(s)
 * and the relative 2-norm of this vector. */
void primalResidual_s(const double* velocity, NumericsMatrix* H, const double* globalVelocity,
                      const double* w, const double* s, double* out, double* rnorm,
                      const double tol);

/* Same as primalResidual_s with a chosen norm type. */
void primalResidual_s_type(const double* velocity, NumericsMatrix* H,
                           const double* globalVelocity, const double* w, const double* s,
                           double* out, double* rnorm, const double tol, const int type);

#ifdef __cplusplus
}
#endif

#endif
