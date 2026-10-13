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

#ifndef fc2d_compute_error_H
#define fc2d_compute_error_H

/*!\file fc2d_compute_error.h
  \brief functions related to error computation for friction-contact problems

*/

#include "NumericsFwd.h"    // for FrictionContactProblem
#include "SolverOptions.h"  // for SolverOptions
#if defined(__cplusplus)
extern "C" {
#endif

void fc2d_unitary_compute_and_add_error(const double* restrict r, const double* restrict u,
                                        double mu, double* restrict error, double* worktmp);

/** Build the modified local velocity used by the 2D natural map.

    For each contact i:
    \f$ u^{mod}_{N,i} = u_{N,i} + \mu_i |u_{T,i}| \f$
    \f$ u^{mod}_{T,i} = u_{T,i} \f$

    \param[in]  nc          number of contacts
    \param[in]  u           velocity vector of length 2*nc
    \param[in]  mu          friction coefficients, array of length nc
    \param[out] u_modified  modified velocity vector of length 2*nc
    \param[in]  work        scratch array of length nc (may be NULL)
 */
void fc2d_compute_modified_local_velocity(int nc, const double* u, const double* mu,
                                          double* u_modified, double* work);

/** Compute the 2D complementarity error with selectable scaling.

    The function mirrors fc3d_compute_complementarity_error for the 2D case.

    \param[in] nc         number of contacts
    \param[in] r          reaction vector of length 2*nc
    \param[in] u          velocity vector of length 2*nc
    \param[in] mu         friction coefficients, array of length nc
    \param[in] norm       scaling norm (see SICONOS_FRICTION_2D_SCALING)
    \param[in] scope      scaling scope (see SICONOS_FRICTION_2D_SCALING_SCOPE)
    \param[in] reduction  reduction type (see SICONOS_FRICTION_2D_SCALING_REDUCTION)
    \param[in] alpha      weight for the weighted cone norm
    \param[in] work       scratch array of length nc
    \param[in,out] tmp    working vector of length 2*nc
    \return the scaled complementarity error
 */
double fc2d_compute_complementarity_error(int nc, const double* r, const double* u,
                                          const double* mu, int norm, int scope, int reduction,
                                          double alpha, double* work, double* tmp);

/** Error computation for friction-contact 2D problem
    \param problem the structure which defines the friction-contact problem
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options solver options (scaling options are read from iparam/dparam)
    \param[in,out] error value
    \return 0 if ok
 */
int fc2d_compute_error(FrictionContactProblem* problem, double* z, double* w, double tolerance,
                       SolverOptions* options, double* error);

#if defined(__cplusplus)
}
#endif

#endif
