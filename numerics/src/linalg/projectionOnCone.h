/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef ProjectionOnCone_H
#define ProjectionOnCone_H

/*!\file projectionOnCone.h
 * \brief function to project on cones
 */

enum { PROJCONE_DUAL, PROJCONE_INSIDE, PROJCONE_BOUNDARY };

#if defined(__cplusplus)
extern "C" {
#endif

/**
   Projection on the second Order Cone in \f$ R^3 \f$

   \f$ K \{r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1  \} \f$

   \param[in,out] r the vector to be projected
   \param[in] mu the angle of the cone
   \return the type of projection
*/
unsigned projectionOnCone(double* r, double mu);

/**
   Natural map associated with the projection on the second order cone.

   Computes \f$ F(r,u) = r - P_K(r - \rho u) \f$ where \f$ P_K \f$ is the
   projection onto the cone \f$ K = \{ r \mid r_1 \geq 0, \sqrt(r_2^2+r_3^2)
   \leq \mu r_1 \} \f$.

   \param[in] r the reaction vector of size 3
   \param[in] u the velocity vector of size 3
   \param[in] rho the scaling parameter (must be positive)
   \param[in] mu the friction coefficient (angle of the cone)
   \param[out] F the natural-map residual of size 3
*/
void naturalMapOnCone(const double r[3], const double u[3], double rho, double mu, double F[3]);

/**
   Compute an element of the the subdifferential of the
   projection on the second Order Cone in \f$ R^3 \f$

   \f$ K \{ r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq mu r_1  \} \f$

   \param[out] H an element of the the subdifferential
   \param[in] r the vector to be projected
   \param[in] mu the angle of the cone
   \return the type of projection
*/
unsigned subdifferentialProjectionOnCone(double* H, double* r, double mu);

/**
    Projection on the second Order Cone in \f$ R^3 \f$

    \f$ K \{ r, r_1 \geq 0, 0 mu \sqrt(u_2^2+u_3^2) \geq u_1  \} \f$

    \param[in,out] u the vector to be projected
    \param[in] mu the angle of the cone
    \return the type of projection
*/
unsigned projectionOnDualCone(double* u, double mu);

/**
   Projection on the second Order Cone in \f$ R^n \f$

   \f$ K \{r, r_1 \geq 0, 0 \|[r_2,r_n]\| \geq mu r_1  \} \f$

   \param[in,out] r the vector to be projected
   \param[in] mu the angle of the cone
   \param[in] size dimension
*/
void projectionOnSecondOrderCone(double* r, double mu, int size);

/**
   Vectorized projection on a stack of 3D friction cones.

   Projects each 3-block of \p r onto the cone \f$ K_i = \
   \{ (r_N,r_T) \mid \|r_T\| \le \mu_i r_N \} \f$ in place.

   \param[in]     nc   number of cones / contacts
   \param[in,out] r    vector of length 3*nc, replaced by its projection
   \param[in]     mu   friction coefficients, array of length nc
   \param[in,out] work scratch array of length nc (may be NULL)
*/
void projectionOnConeVector(int nc, double* r, const double* mu, double* work);

/**
   Vectorized natural map for a stack of 3D friction cones.

   Computes \f$ F_i(r_i,u_i) = r_i - P_{K_i}(r_i - \rho u_i) \f$ for
   \f$ i = 0 \ldots nc-1 \f$, where each block has size 3 and \f$ K_i \f$ is the
   cone with friction coefficient mu[i].

   A workspace of length nc must be provided in work. If work is NULL, the
   function allocates and frees a temporary buffer internally (slower if called
   repeatedly).

   \param[in] nc number of cones / contacts
   \param[in] r reaction vector of length 3*nc
   \param[in] u velocity vector of length 3*nc
   \param[in] rho scaling parameter (must be positive)
   \param[in] mu friction coefficients, array of length nc
   \param[out] F natural-map residual vector of length 3*nc
   \param[in,out] work scratch array of length nc (may be NULL)
*/
void naturalMapOnConeVector(int nc, const double* r, const double* u, double rho,
                            const double* mu, double* F, double* work);

#if defined(__cplusplus)
}
#endif

#endif
