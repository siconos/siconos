// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef projection_soc_H
#define projection_soc_H

/*!\file projection_soc.h
  \brief Euclidean projection on second-order (friction) cones and its Jacobian.

  The cone of coefficient \f$ \mu \geq 0 \f$ in \f$ \mathbb{R}^n \f$ is
  \f[
    K_\mu = \{ x = (x_N, x_T) \in \mathbb{R} \times \mathbb{R}^{n-1}
               \mid \|x_T\| \leq \mu x_N \}
  \f]
  (Coulomb friction cone, Drucker-Prager cone with \f$ \mu = \eta \f$). Its
  dual cone is \f$ K_\mu^* = K_{1/\mu} \f$. With \f$ t = \|x_T\| \f$ the projection is
  - \f$ x \f$ if \f$ t \leq \mu x_N \f$ (inside, ::SOC_PROJECTION_INSIDE);
  - \f$ 0 \f$ if \f$ \mu t \leq -x_N \f$ (polar cone, ::SOC_PROJECTION_POLAR);
  - \f$ \alpha (1, \mu x_T / t) \f$, \f$ \alpha = (x_N + \mu t)/(1 + \mu^2) \f$
    otherwise (boundary, ::SOC_PROJECTION_BOUNDARY).

  Natural maps and error computations built on these projections are in
  convex_set.h and natural_map_error.h.
*/

/** Region of the space where the projected point was (return value of the projections) */
enum SOC_PROJECTION_STATUS {
  /** the point is in the polar cone: projected on the apex 0 */
  SOC_PROJECTION_POLAR = 0,
  /** the point is in the cone: unchanged */
  SOC_PROJECTION_INSIDE = 1,
  /** the point is projected on the boundary of the cone */
  SOC_PROJECTION_BOUNDARY = 2
};

#if defined(__cplusplus)
extern "C" {
#endif

/** Projection on the cone \f$ K_\mu \f$ of \f$ \mathbb{R}^3 \f$.

    \param[in,out] x the vector of size 3 to be projected
    \param[in] mu the coefficient of the cone
    \return the region (see ::SOC_PROJECTION_STATUS)
*/
unsigned soc_projection_3d(double* x, double mu);

/** Projection on the cone \f$ K_\mu \f$ of \f$ \mathbb{R}^2 \f$
    (\f$ |x_T| \leq \mu x_N \f$).

    \param[in,out] x the vector of size 2 to be projected
    \param[in] mu the coefficient of the cone
    \return the region (see ::SOC_PROJECTION_STATUS)
*/
unsigned soc_projection_2d(double* x, double mu);

/** Projection on the cone \f$ K_\mu \f$ of \f$ \mathbb{R}^n \f$, any \f$ n \geq 1 \f$
    (\f$ n = 1 \f$: projection on \f$ \mathbb{R}_+ \f$).

    \param[in,out] x the vector of size dim to be projected
    \param[in] mu the coefficient of the cone
    \param[in] dim the dimension n
    \return the region (see ::SOC_PROJECTION_STATUS)
*/
unsigned soc_projection(double* x, double mu, int dim);

/** Projection on the dual cone \f$ K_\mu^* = \{ x \mid \mu \|x_T\| \leq x_N \} \f$
    of \f$ \mathbb{R}^3 \f$.

    \param[in,out] x the vector of size 3 to be projected
    \param[in] mu the coefficient of the (primal) cone
    \return the region (see ::SOC_PROJECTION_STATUS)
*/
unsigned soc_dual_projection_3d(double* x, double mu);

/** An element of the (Clarke) generalized Jacobian of the projection on
    \f$ K_\mu \f$ in \f$ \mathbb{R}^3 \f$.

    \param[out] J the 3x3 Jacobian, column-major
    \param[in] x the point where the projection is differentiated (not modified)
    \param[in] mu the coefficient of the cone
    \return the region (see ::SOC_PROJECTION_STATUS)
*/
unsigned soc_projection_jacobian_3d(double* J, const double* x, double mu);

/** Projection on a stack of cones of \f$ \mathbb{R}^3 \f$, in place:
    \f$ x_i \leftarrow P_{K_{\mu_i}}(x_i) \f$, \f$ i = 0 \ldots n_c-1 \f$.

    \param[in] nc number of cones
    \param[in,out] x vector of length 3 nc
    \param[in] mu coefficients of the cones, length nc
    \param[in,out] work scratch array of length nc (may be NULL: allocated internally)
*/
void soc_projection_3d_stack(int nc, double* x, const double* mu, double* work);

#if defined(__cplusplus)
}
#endif

#endif
