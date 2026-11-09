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

#ifndef fc2d_AlartCurnier_functions_h
#define fc2d_AlartCurnier_functions_h

/*!\file fc2d_AlartCurnier_functions.h
 * \brief Alart--Curnier function and Jacobians for a single 2D frictional
 *        contact.
 *
 * For a 2D contact the unknowns are the local velocity
 * \f$u = (u_N, u_T)\f$ and the reaction \f$r = (r_N, r_T)\f$.  The
 * Alart--Curnier residual is
 *
 * \f[
 *   \Phi_{\mathrm{AC}}(u, r, \rho) =
 *   \begin{pmatrix}
 *     \rho_N u_N & \text{if } r_N - \rho_N u_N \ge 0, \\
 *     r_N        & \text{otherwise},
 *   \end{pmatrix}
 * \f]
 * for the normal part and
 *
 * \f[
 *   \Phi_{\mathrm{AC},T}(u, r, \rho) =
 *   \begin{cases}
 *     \rho_T u_T & \text{if } |r_T - \rho_T u_T| \le \mu r_N, \\
 *     r_T - \mu r_N \operatorname{sgn}(r_T - \rho_T u_T) & \text{otherwise},
 *   \end{cases}
 * \f]
 * for the tangential part (Jean--Moreau radius \f$\mu r_N\f$).
 */

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * \brief Compute the 2D Alart--Curnier function and its Jacobians.
 *
 * \param[in] reaction contact reaction \f$r = (r_N, r_T)\f$
 * \param[in] velocity local velocity \f$u = (u_N, u_T)\f$
 * \param[in] mu friction coefficient
 * \param[in] rho scaling parameters \f$\rho = (\rho_N, \rho_T)\f$
 * \param[out] residual Alart--Curnier residual \f$\Phi_{\mathrm{AC}}\f$
 *                 (size 2); may be NULL
 * \param[out] jacobian_u partial Jacobian \f$\partial \Phi / \partial u\f$
 *                 stored in column-major format (size 2x2); may be NULL
 * \param[out] jacobian_r partial Jacobian \f$\partial \Phi / \partial r\f$
 *                 stored in column-major format (size 2x2); may be NULL
 */
void computeAlartCurnier2D(double* reaction, double* velocity, double mu, double* rho,
                           double* residual, double* jacobian_u, double* jacobian_r);
void computeAlartCurnier2D_reformulated(double* reaction, double* velocity, double mu,
                                        double* rho, double* residual, double* jacobian_u,
                                        double* jacobian_r);
void computeAlartCurnier2D_AC(double* reaction, double* velocity, double mu, double* rho,
                              double* residual, double* jacobian_u, double* jacobian_r);
void computeAlartCurnier2D_JeanMoreau(double* reaction, double* velocity, double mu,
                                      double* rho, double* residual, double* jacobian_u,
                                      double* jacobian_r);

#if defined(__cplusplus)
}
#endif

#endif
