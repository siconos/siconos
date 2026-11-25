// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NONSMOOTH_NEWTON_LEVENBERG_MARQUARDT_H
#define NONSMOOTH_NEWTON_LEVENBERG_MARQUARDT_H

/*!\file nonsmooth_newton_levenberg_marquardt.h
 * \brief Levenberg--Marquardt fallback of the semi-smooth Newton solvers.
 *
 * The semi-smooth Newton solvers of the friction problems (GFC3D_NSN_AC, GFC2D_NSN_AC,
 * FC3D_NSN_AC/FB/NM) solve \f$ J(x_k) d_k = -\Psi(x_k) \f$ and globalize the iteration with a
 * line search on the merit function \f$ q(x) = \frac{1}{2}\|\Psi(x)\|^2 \f$. When the Jacobian
 * is (nearly) singular, the Newton direction (or its shifted variant
 * \f$ (J + \epsilon I) d = -\Psi \f$) has huge components along the directions of the small
 * singular values of J; the line search then only accepts tiny steps and the iteration stalls
 * at a point where \f$ \Psi \neq 0 \f$.
 *
 * The fallback is not a universal improvement: on the 2D granular problems of the GFC2D test
 * collection, a full Levenberg--Marquardt step moved the iterate to a wrong region (almost all
 * the contacts switched to sliding) from which the iteration did not recover, while Newton
 * alone converged after a slow phase. It is therefore enabled by default only for GFC3D_NSN_AC
 * (SICONOS_FRICTION_3D_NSN_LM_FALLBACK = 1); it is available, disabled by default, for
 * GFC2D_NSN_AC (SICONOS_FRICTION_2D_NSN_LM_FALLBACK) and FC3D_NSN_AC/FB/NM.
 *
 * On a stalled iteration (see nonsmooth_newton_lm_is_stalled) the solvers use instead the
 * Levenberg--Marquardt direction
 * \f[ (J^\top J + \lambda I) d = -J^\top \Psi, \qquad \lambda = c \|\Psi\|, \f]
 * which is a descent direction of q for any \f$ \lambda > 0 \f$ and damps the components along
 * the small singular values \f$ \sigma \f$ of J (factor \f$ \sigma / (\sigma^2 + \lambda) \f$)
 * instead of amplifying them. Since \f$ \lambda \to 0 \f$ with \f$ \|\Psi\| \f$, the fast
 * local convergence of Newton is recovered as soon as the iteration progresses again.
 *
 * See docs/devel/GFC3D_NSN_LM_FALLBACK.md.
 */

#include "NumericsFwd.h"  // for NumericsMatrix

#if defined(__cplusplus)
extern "C" {
#endif

/** Relative decrease of the residual norm below which a Newton iteration does not progress:
 *  \f$ \|\Psi_k\| > \f$ NONSMOOTH_NEWTON_LM_STALL_RATIO \f$ \|\Psi_{k-1}\| \f$. */
#define NONSMOOTH_NEWTON_LM_STALL_RATIO 0.99

/** Number of consecutive iterations without progress after which the Newton iteration is
 *  considered stalled. With 1 (validated on the GFC3D FCLIB problems) the fallback rescues the
 *  stalls of GFC3D_NSN_AC; with 3 it came too late on these problems. */
#define NONSMOOTH_NEWTON_LM_STALL_ITERATIONS 1

/** Stall test of the Levenberg--Marquardt fallback.
 *  \param residual_norm \f$ \|\Psi_k\| \f$
 *  \param previous_residual_norm \f$ \|\Psi_{k-1}\| \f$ (INFINITY at the first iteration)
 *  \param[in,out] stalled_iterations number of consecutive iterations without progress
 *  (initialized to 0 by the caller, updated here)
 *  \return 1 if the residual decreased by less than 1% at each of the last
 *  NONSMOOTH_NEWTON_LM_STALL_ITERATIONS iterations
 */
int nonsmooth_newton_lm_is_stalled(double residual_norm, double previous_residual_norm,
                                   int* stalled_iterations);

/** Levenberg--Marquardt direction: solves \f$ (J^\top J + \lambda I) d = -J^\top \Psi \f$.
 *
 * J may be stored as NM_DENSE, NM_SPARSE_BLOCK or NM_SPARSE (a sparse copy is used for the
 * dense and sparse-block storages). The normal matrix is factorized with an LU solve.
 *
 * \param[in] jacobian the (square) Jacobian J
 * \param[in] residual the residual \f$ \Psi \f$
 * \param[in] lambda the damping parameter \f$ \lambda > 0 \f$
 * \param[out] direction the direction d (size of J)
 * \return the info of the linear solver (0 on success)
 */
int nonsmooth_newton_lm_direction(NumericsMatrix* jacobian, const double* residual,
                                  double lambda, double* direction);

#if defined(__cplusplus)
}
#endif

#endif
