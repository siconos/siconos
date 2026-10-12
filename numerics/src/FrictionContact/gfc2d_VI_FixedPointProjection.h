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

#ifndef GFC2D_VI_FIXEDPOINTPROJECTION_H
#define GFC2D_VI_FIXEDPOINTPROJECTION_H

/*!\file gfc2d_VI_FixedPointProjection.h
 * \brief VI fixed point projection solver for 2D global friction contact problems.
 */

#include "NumericsFwd.h"  // for GlobalFrictionContactProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * \brief VI Fixed Point Projection solver for 2D global friction contact problems.
 *
 * The solver follows the same strategy as gfc3d_VI_FixedPointProjection:
 * the global 2D friction-contact problem is reformulated as a variational
 * inequality on the mixed variable (globalVelocity, reaction) and solved
 * with a fixed-point projection scheme.  The reaction part is projected
 * onto the 2D Coulomb cone.
 *
 * \param[in] problem the global 2D friction contact problem
 * \param[in,out] reaction contact reaction (size: 2 * numberOfContacts)
 * \param[out] velocity local contact velocity (size: 2 * numberOfContacts)
 * \param[in,out] globalVelocity global velocity (size: M->size0)
 * \param[out] info convergence info (0 if converged)
 * \param[in,out] options solver options
 */
void gfc2d_VI_FixedPointProjection(GlobalFrictionContactProblem *problem, double *reaction,
                                   double *velocity, double *globalVelocity, int *info,
                                   SolverOptions *options);

#if defined(__cplusplus)
}
#endif

#endif
