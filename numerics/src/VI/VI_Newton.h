// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef VI_NEWTON_H
#define VI_NEWTON_H

/*!\file VI_Newton.h
 * \brief Functions for solving VI using Newton method
 */

#include "NumericsMatrix.h"  // IWYU pragma: keep

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
extern "C" {
#endif

void VI_compute_F(void* data_opaque, double* x, double* F);
void VI_compute_error_box(void* data_opaque, double* x, double* F, double* Jac_F_merit,
                          double tol, double* err);
void VI_compute_F_box_Qi(void* data_opaque, double* x, double* F, double* Fbox);
void VI_compute_H_box_Qi(void* data_opaque, double* x, double* F, double* workV1,
                         double* workV2, NumericsMatrix* H);

#if defined(__cplusplus)
}
#endif

#endif
