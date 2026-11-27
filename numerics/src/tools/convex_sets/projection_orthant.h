// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef projection_orthant_H
#define projection_orthant_H

/*!\file projection_orthant.h
  \brief Euclidean projection on the nonnegative orthant and on boxes, and its Jacobian.

  - orthant \f$ \mathbb{R}^n_+ \f$ (LCP, MLCP complementarity rows):
    \f$ P(x)_i = \max(0, x_i) \f$;
  - box \f$ [l, u] = \{ x \mid l_i \leq x_i \leq u_i \} \f$ (Relay, bounds):
    \f$ P(x)_i = \min(u_i, \max(l_i, x_i)) \f$.

  The Jacobian is diagonal with entries 0 or 1: only its diagonal is returned.
*/

#if defined(__cplusplus)
extern "C" {
#endif

/** Projection on \f$ \mathbb{R}^n_+ \f$, in place.

    \param[in,out] x vector of size n
    \param[in] n dimension
*/
void orthant_projection(double* x, int n);

/** Diagonal of an element of the generalized Jacobian of the projection on
    \f$ \mathbb{R}^n_+ \f$: 1 if \f$ x_i > 0 \f$, 0 otherwise.

    \param[out] d diagonal of the Jacobian, size n
    \param[in] x the point
    \param[in] n dimension
*/
void orthant_projection_jacobian_diag(double* d, const double* x, int n);

/** Projection on the box \f$ [lb, ub] \f$, in place.

    \param[in,out] x vector of size n
    \param[in] n dimension
    \param[in] lb lower bounds, size n
    \param[in] ub upper bounds, size n
*/
void box_projection(double* x, int n, const double* lb, const double* ub);

/** Diagonal of an element of the generalized Jacobian of the projection on
    \f$ [lb, ub] \f$: 1 if \f$ lb_i < x_i < ub_i \f$, 0 otherwise.

    \param[out] d diagonal of the Jacobian, size n
    \param[in] x the point
    \param[in] n dimension
    \param[in] lb lower bounds, size n
    \param[in] ub upper bounds, size n
*/
void box_projection_jacobian_diag(double* d, const double* x, int n, const double* lb,
                                  const double* ub);

#if defined(__cplusplus)
}
#endif

#endif
