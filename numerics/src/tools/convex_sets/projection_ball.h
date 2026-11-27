// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef projection_ball_H
#define projection_ball_H

/*!\file projection_ball.h
  \brief Euclidean projection on balls and cylinders, and its Jacobian.

  Sets of \f$ \mathbb{R}^n \f$, \f$ R \geq 0 \f$, \f$ x = (x_N, x_T) \in \mathbb{R}
  \times \mathbb{R}^{n-1} \f$:
  - ball \f$ B_R = \{ x \mid \|x\| \leq R \} \f$ (disk for \f$ n = 2 \f$, Von Mises
    plane-stress criterion in scaled variables);
  - cylinder \f$ C_R = \mathbb{R} \times B_R^{n-1} = \{ x \mid \|x_T\| \leq R \} \f$
    (Von Mises plane-strain criterion: no condition on the hydrostatic part);
  - positive cylinder \f$ C_R^+ = \mathbb{R}_+ \times B_R^{n-1} \f$ (Tresca friction).

  The projection on a ball is the radial return \f$ R x / \|x\| \f$ outside the ball;
  the projection on a product set is the product of the projections.
*/

/** Region of the space where the projected point was (return value of the projections) */
enum BALL_PROJECTION_STATUS {
  /** the point is in the set: unchanged */
  BALL_PROJECTION_INSIDE = 1,
  /** the point is projected on the boundary of the set */
  BALL_PROJECTION_BOUNDARY = 2
};

#if defined(__cplusplus)
extern "C" {
#endif

/** Projection on the ball \f$ B_R \f$ of \f$ \mathbb{R}^{dim} \f$.

    \param[in,out] x vector of size dim
    \param[in] dim dimension
    \param[in] R radius
    \return the region (see ::BALL_PROJECTION_STATUS)
*/
unsigned ball_projection(double* x, int dim, double R);

/** An element of the generalized Jacobian of the projection on \f$ B_R \f$:
    \f$ I \f$ inside, \f$ (R/\|x\|)(I - x x^T/\|x\|^2) \f$ outside.

    \param[out] J dim x dim Jacobian, column-major
    \param[in] x the point (not modified)
    \param[in] dim dimension
    \param[in] R radius
    \return the region (see ::BALL_PROJECTION_STATUS)
*/
unsigned ball_projection_jacobian(double* J, const double* x, int dim, double R);

/** Projection on the cylinder \f$ C_R = \mathbb{R} \times B_R^{dim-1} \f$: \f$ x_N \f$
    unchanged, \f$ x_T \f$ projected on the ball.

    \param[in,out] x vector of size dim
    \param[in] dim dimension
    \param[in] R radius
    \return the region (see ::BALL_PROJECTION_STATUS)
*/
unsigned cylinder_projection(double* x, int dim, double R);

/** An element of the generalized Jacobian of the projection on \f$ C_R \f$.

    \param[out] J dim x dim Jacobian, column-major
    \param[in] x the point (not modified)
    \param[in] dim dimension
    \param[in] R radius
    \return the region (see ::BALL_PROJECTION_STATUS)
*/
unsigned cylinder_projection_jacobian(double* J, const double* x, int dim, double R);

/** Projection on the positive cylinder \f$ C_R^+ = \mathbb{R}_+ \times B_R^{dim-1} \f$
    (Tresca friction): \f$ x_N \leftarrow \max(0, x_N) \f$, \f$ x_T \f$ projected on the
    ball.

    \param[in,out] x vector of size dim
    \param[in] dim dimension
    \param[in] R radius
    \return the region (see ::BALL_PROJECTION_STATUS)
*/
unsigned positive_cylinder_projection(double* x, int dim, double R);

#if defined(__cplusplus)
}
#endif

#endif
