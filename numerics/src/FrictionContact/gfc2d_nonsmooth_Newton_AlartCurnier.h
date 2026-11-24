/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef gfc2d_nonsmooth_Newton_AlartCurnier_h
#define gfc2d_nonsmooth_Newton_AlartCurnier_h

/*!\file gfc2d_nonsmooth_Newton_AlartCurnier.h
 *
 * \brief Public interface of the sparse global Alart--Curnier non-smooth
 *        Newton solver for 2D frictional contact problems.
 *
 * The solver addresses the global 2D frictional contact problem
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
 */

#include "NumericsFwd.h"  // for SolverOptions, GlobalFrictionContactProblem

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * \brief Solve a sparse global 2D friction contact problem with the
 *        Alart--Curnier non-smooth Newton method.
 *
 * \param[in] problem the global friction contact problem
 * \param[in,out] reaction on input, initial guess for the contact reaction;
 *                 on output, the computed reaction (size: 2 x numberOfContacts)
 * \param[in,out] velocity on input, initial guess for the local contact
 *                 velocity; on output, the computed velocity
 *                 (size: 2 x numberOfContacts)
 * \param[in,out] globalVelocity on input, initial guess for the global
 *                 velocity; on output, the computed global velocity
 *                 (size: M->size0)
 * \param[out] info solver status: 0 on success, nonzero otherwise
 * \param[in,out] options solver options and parameters
 */
void gfc2d_nonsmooth_Newton_AlartCurnier(GlobalFrictionContactProblem* problem,
                                         double* reaction, double* velocity,
                                         double* globalVelocity, int* info,
                                         SolverOptions* options);

/**
 * \brief Run a VI-FPP warm start for GFC2D_NSN_AC.
 *
 * This function can be called before the Newton solver (it is called
 * automatically by gfc2d_nonsmooth_Newton_AlartCurnier when the option
 * SICONOS_FRICTION_2D_NSN_WARM_START_STRATEGY is set to YES) or as a
 * standalone fallback/preconditioner.
 *
 * \param[in] problem the global 2D friction contact problem
 * \param[in,out] reaction initial guess / warm-started reaction
 * \param[in,out] velocity initial guess / warm-started local velocity
 * \param[in,out] globalVelocity initial guess / warm-started global velocity
 * \param[in,out] options GFC2D_NSN_AC solver options; on output the rho
 *                        parameter may be overwritten with the VI-FPP rho
 * \return 0 if the warm start ran (or was disabled), nonzero on error
 */
int gfc2d_nsn_ac_warm_start(GlobalFrictionContactProblem* problem, double* reaction,
                            double* velocity, double* globalVelocity, SolverOptions* options);

/**
 * \brief Set default solver options for GFC2D_NSN_AC.
 *
 * \param[in,out] options the solver options structure to initialize
 */
void gfc2d_nsn_ac_set_default(SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
