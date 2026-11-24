// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MIN_MERIT
#define MIN_MERIT
/*!\file min_merit.h

  \brief Functions for the min based merit function

  A set of routines used in the min reformulation of a CP

  The min function is:
  \f{equation*}
  \mathbf{F}_{\mathrm{min}}(z) = \min( z, F(z))
  \f}
*/

#include "NumericsFwd.h"  // for NumericsMatrix

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/** Compute \f$ \mathbf{F}_{\mathrm{min}}(z) \f$, possibly in the mixed case
 * \param n1 number of equality constraints
 * \param n2 size of the complementary variables
 * \param[in] z input vector
 * \param[in] F value of F
 * \param[out] Fmin returned vector
 */
void F_min(int n1, int n2, double* restrict z, double* restrict F, double* restrict Fmin);

/** Compute an element of Jac F_min
 * \param n1 number of equality constraints
 * \param n2 size of the complementarity variables
 * \param[in] z input vector
 * \param[in] F value of F
 * \param[in] nabla_F value of nabla_F
 * \param[out] H returned vector
 */
void Jac_F_min(int n1, int n2, double* restrict z, double* restrict F, NumericsMatrix* nabla_F,
               NumericsMatrix* H);

#if defined(__cplusplus)
}
#endif

#endif
