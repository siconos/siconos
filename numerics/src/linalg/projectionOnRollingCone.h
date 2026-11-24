// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef ProjectionOnRollingCone_H
#define ProjectionOnRollingCone_H

/*!\file projectionOnCone.h
 * function to project on rolling cones
 */

enum {
  PROJRCONE_DUAL,
  PROJRCONE_INSIDE,
  PROJRCONE_BOUNDARY_FRICTION,
  PROJRCONE_BOUNDARY_ROLLING,
  PROJRCONE_BOUNDARY_FRICTION_ROLLING
};
#if defined(__cplusplus)
extern "C" {
#endif

/**
    projectionOnRollingCone Projection on the second Order Cone in \f$ R^3 \f$, \f$ K
   \{ r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \}
   \f$

    \param[in,out] r the vector to be projected
    \param[in] mu the angle of the cone in reaction
    \param[in] mu_r the angle of the cone in moment
    \return the type of projection
*/
unsigned int projectionOnRollingCone(double* r, double mu, double mur);

/**
   projectionOnRollingCone Projection on the second Order Cone in \f$ R^3 \f$, \f$ K \{
   r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \} \f$

   \param[in,out] r the vector to be projected
   \param[in] mu the angle of the cone in reaction
   \param[in] mu_r the angle of the cone in moment
   \return the type of projection
*/
unsigned int projectionOn2DRollingCone(double* r, double mu, double mur);

/**
   projectionOnDualRollingCone Projection on the second Order Cone in \f$ R^3 \f$, \f$ K
   \{ r, r_1 \geq 0, 0 mu \sqrt(u_2^2+u_3^2) \geq u_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \}
   \f$

   \param[in,out] u the vector to be projected
   \param[in] mu the angle of the cone in reaction
   \param[in] mu_r the angle of the cone in moment
   \return the type of projection
*/
unsigned projectionOnDualRollingCone(double* u, double mu, double mur);

void display_status_rolling_cone(unsigned int status);

/**
   subdifferentialProjectionOnRollingCone.
   Compute an element of the the subdifferential of the
   projection on the second Order Cone in \f$ R^5 \f$, \f$ K \{
   r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1, 0 \sqrt(r_3^2+r_4^2) \geq mu_r r_1  \} \f$

   \param[out] H an element of the the subdifferential
   \param[in] r the vector to be projected
   \param[in] mu the angle of the cone in reaction
   \param[in] mu_r the angle of the cone in moment
   \return the type of projection
*/

unsigned subdifferentialProjectionOnRollingCone(double* H, double* r, double mu, double mur);

#if defined(__cplusplus)
}
#endif

#endif
