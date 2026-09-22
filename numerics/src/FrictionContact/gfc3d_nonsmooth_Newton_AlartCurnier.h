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

#ifndef gfc3d_nonsmooth_Newton_AlartCurnier_h
#define gfc3d_nonsmooth_Newton_AlartCurnier_h

/*!\file gfc3d_nonsmooth_Newton_AlartCurnier.h
 *
 * \brief Public interface of the sparse global Alart--Curnier non-smooth
 *        Newton solver for 3D frictional contact problems.
 *
 * The solver addresses the global 3D frictional contact problem
 *
 * \f[
 *   \begin{cases}
 *     M \, v = H \, r + q, \\
 *     u = H^\top v + b, \\
 *     (u_N, r_N, r_T) \text{ satisfy Signorini--Coulomb conditions},
 *   \end{cases}
 * \f]
 *
 * where \f$v\f$ is the global velocity, \f$u\f$ the local contact velocity,
 * \f$r\f$ the contact reaction, \f$M\f$ the mass matrix, \f$H\f$ the
 * contact Jacobian and \f$q, b\f$ data vectors.
 *
 * The contact conditions are reformulated as a semi-smooth residual using the
 * Alart--Curnier function.  A non-smooth Newton iteration is applied to the
 * concatenated residual
 *
 * \f[
 *   \Psi(v, u, r) =
 *   \begin{pmatrix}
 *     M v - H r - q \\
 *     u - H^\top v - b \\
 *     \Phi_{\mathrm{AC}}(u, r, \rho)
 *   \end{pmatrix},
 * \f]
 *
 * with \f$\Phi_{\mathrm{AC}}\f$ the Alart--Curnier residual and
 * \f$\rho\f$ a vector of contact scaling parameters.
 */

#include "NumericsFwd.h"  // for SolverOptions, GlobalFrictionContactProblem

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * \brief Solve a sparse global 3D friction contact problem with the
 *        Alart--Curnier non-smooth Newton method.
 *
 * \param[in] problem the global friction contact problem
 * \param[in,out] reaction on input, initial guess for the contact reaction;
 *                 on output, the computed reaction (size: 3 x numberOfContacts)
 * \param[in,out] velocity on input, initial guess for the local contact
 *                 velocity; on output, the computed velocity
 *                 (size: 3 x numberOfContacts)
 * \param[in,out] globalVelocity on input, initial guess for the global
 *                 velocity; on output, the computed global velocity
 *                 (size: M->size0)
 * \param[out] info solver status: 0 on success, nonzero otherwise
 * \param[in,out] options solver options and parameters
 */
void gfc3d_nonsmooth_Newton_AlartCurnier(GlobalFrictionContactProblem* problem,
                                         double* reaction, double* velocity,
                                         double* globalVelocity, int* info,
                                         SolverOptions* options);

/**
 * \brief Set default solver options for GFC3D_NSN_AC.
 *
 * \param[in,out] options the solver options structure to initialize
 */
void gfc3d_nsn_ac_set_default(SolverOptions* options);

/**
 * \brief No-op initialization hook for the sparse global Alart--Curnier
 *        solver.
 *
 * This function is part of the solver registration interface.  It currently
 * performs no operation because all internal allocations are handled lazily
 * inside gfc3d_nonsmooth_Newton_AlartCurnier.
 *
 * \param[in,out] options the solver options structure
 */
void gfc3d_sparseGlobalAlartCurnierInit(SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
