// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef projection_rolling_cone_H
#define projection_rolling_cone_H

/*!\file projection_rolling_cone.h
  \brief Euclidean projection on rolling friction cones and its Jacobian.

  Rolling friction cone in \f$ \mathbb{R}^5 \f$ (3D: normal, 2 tangential, 2 rolling
  components) or \f$ \mathbb{R}^3 \f$ (2D: normal, tangential, rolling):
  \f$ K = \{ r \mid \|r_T\| \leq \mu r_N, \|r_R\| \leq \mu_r r_N \} \f$.
*/

/** Region of the space where the projected point was (return value of the projections) */
enum ROLLING_CONE_PROJECTION_STATUS {
  ROLLING_CONE_PROJECTION_POLAR,
  ROLLING_CONE_PROJECTION_INSIDE,
  ROLLING_CONE_PROJECTION_BOUNDARY_FRICTION,
  ROLLING_CONE_PROJECTION_BOUNDARY_ROLLING,
  ROLLING_CONE_PROJECTION_BOUNDARY_FRICTION_ROLLING
};
#if defined(__cplusplus)
extern "C" {
#endif

/**
    rolling_cone_projection Projection on the second Order Cone in \f$ R^3 \f$, \f$ K
   \{ r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \}
   \f$

    \param[in,out] r the vector to be projected
    \param[in] mu the angle of the cone in reaction
    \param[in] mu_r the angle of the cone in moment
    \return the type of projection
*/
unsigned int rolling_cone_projection(double* r, double mu, double mur);

/**
   rolling_cone_projection Projection on the second Order Cone in \f$ R^3 \f$, \f$ K \{
   r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \} \f$

   \param[in,out] r the vector to be projected
   \param[in] mu the angle of the cone in reaction
   \param[in] mu_r the angle of the cone in moment
   \return the type of projection
*/
unsigned int rolling_cone_2d_projection(double* r, double mu, double mur);

/**
   rolling_cone_dual_projection Projection on the second Order Cone in \f$ R^3 \f$, \f$ K
   \{ r, r_1 \geq 0, 0 mu \sqrt(u_2^2+u_3^2) \geq u_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \}
   \f$

   \param[in,out] u the vector to be projected
   \param[in] mu the angle of the cone in reaction
   \param[in] mu_r the angle of the cone in moment
   \return the type of projection
*/
unsigned rolling_cone_dual_projection(double* u, double mu, double mur);

void rolling_cone_projection_status_display(unsigned int status);

/**
   rolling_cone_projection_jacobian.
   Compute an element of the the subdifferential of the
   projection on the second Order Cone in \f$ R^5 \f$, \f$ K \{
   r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \} \f$

   \param[out] H an element of the the subdifferential
   \param[in] r the vector to be projected
   \param[in] mu the angle of the cone in reaction
   \param[in] mu_r the angle of the cone in moment
   \return the type of projection
*/

unsigned rolling_cone_projection_jacobian(double* H, double* r, double mu, double mur);

#if defined(__cplusplus)
}
#endif

#endif
