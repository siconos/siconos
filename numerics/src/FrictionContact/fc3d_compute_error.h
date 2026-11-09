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

#ifndef fc3d_compute_error_H
#define fc3d_compute_error_H

/*!\file fc3d_compute_error.h
  \brief functions related to error computation for friction-contact problems

*/

#include "FrictionContact_options.h"  // for SICONOS_FRICTION_3D_SCALING
#include "NumericsFwd.h"              // for FrictionContactProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/**
    Error computation using the norm infinity conic for friction-contact 3D problem

    \param problem the structure which defines the friction-contact problem
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param norm norm of a vector (problem->q) for relative error
    \param[in,out] error value
    \return 0 if ok
 */
int fc3d_compute_error_norm_infinity_conic(FrictionContactProblem* problem, double* z,
                                           double* w, double tolerance, SolverOptions* options,
                                           double norm, double* error, int on_dual_cone);

/**
    Error computation on dual cone (using the normal map residual) for one friction-contact 3D
   problem

    \param r the reaction force
    \param u the local velocity
    \param mu coeficient of friction
    \param worktmp work vector
    \param[in,out] error value
 */
void fc3d_unitary_compute_dual_and_add_error(double* restrict r, double* restrict u, double mu,
                                             double* error, double* worktmp);

/**
    Error computation for a friction-contact 3D problem

    \param problem the structure which defines the friction-contact problem
    \param z vector
    \param w vector
    \param options
    \param tolerance value for error computation
    \param[in,out] error value
    \return 0 if ok
 */
int fc3d_compute_error_velocity(FrictionContactProblem* problem, double* z, double* w,
                                double tolerance, SolverOptions* options, double* error);

/**
    Error computation for a friction-contact 3D problem with selectable scaling.

    The scaling is controlled by three options in \p options->iparam:
      - SICONOS_FRICTION_3D_SCALING_NORM     (see SICONOS_FRICTION_3D_SCALING)
      - SICONOS_FRICTION_3D_SCALING_SCOPE    (see SICONOS_FRICTION_3D_SCALING_SCOPE)
      - SICONOS_FRICTION_3D_SCALING_REDUCTION (see SICONOS_FRICTION_3D_SCALING_REDUCTION)
    For SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED the weight is read from
    options->dparam[SICONOS_FRICTION_3D_SCALING_ALPHA].

    \param problem the structure which defines the friction-contact problem
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param[in,out] error value
    \return 0 if ok
*/
int fc3d_compute_error(FrictionContactProblem* problem, double* z, double* w, double tolerance,
                       SolverOptions* options, double* error);

/**
    Compute the modified local velocity for a stack of 3D friction cones.

    For each contact i, the modified velocity is
    \f$ \tilde u_i = ( u_{N,i} + \mu_i \|u_{T,i}\|, u_{T,i} ) \f$.
    The tangential norm is computed with a vectorized loop over all contacts.

    \param[in]  nc          number of contacts
    \param[in]  u           local velocity vector of length 3*nc
    \param[in]  mu          friction coefficients, array of length nc
    \param[out] u_modified  modified local velocity vector of length 3*nc
    \param[in,out] work     scratch array of length nc (may be NULL)
 */
void fc3d_compute_modified_local_velocity(int nc, const double* u, const double* mu,
                                          double* u_modified, double* work);

/**
    Vectorized computation of the (scaled) complementarity error for a stack of 3D
    friction cones.

    Computes the modified local velocity \f$ \tilde u \f$, evaluates the unscaled
    natural map \f$ F(r, \tilde u) = r - P_K(r - \tilde u) \f$, and scales the
    residual according to the selected norm, scope and reduction.
    The function returns \f$ \|F_s\|_2 \f$.

    \param[in] nc          number of contacts
    \param[in] r           reaction vector of length 3*nc
    \param[in] u           local velocity vector of length 3*nc
    \param[in] mu          friction coefficients, array of length nc
    \param[in] norm        norm used for the scale (see SICONOS_FRICTION_3D_SCALING)
    \param[in] scope       per-cone or whole-vector scope
    \param[in] reduction   max, geometric mean or mixed reduction
    \param[in] alpha       weight for SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED
    \param[in,out] work    scratch array of length nc
    \param[in,out] tmp     scratch array of length 3*nc
    \return the scaled complementarity error (Euclidean norm)
 */
double fc3d_compute_complementarity_error(int nc, const double* r, const double* u,
                                          const double* mu, int norm, int scope, int reduction,
                                          double alpha, double* work, double* tmp);

void fc3d_unitary_compute_and_add_error(const double r[3], const double u[3], double mu,
                                        double* error, double worktmp[3]);

/**
    Error computation for one friction-contact 3D problem

    \param z vector
    \param w vector
    \param R radius of the cylinder
    \param worktmp work vector
    \param[in,out] error value
 */
void fc3d_Tresca_unitary_compute_and_add_error(double z[3], double w[3], double R,
                                               double* error, double* worktmp);

/**
    Error computation for friction-contact 3D problem with Tresca Friction

    \param problem the structure which defines the friction-contact problem
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param[in,out] error value
    \return 0 if ok
 */
int fc3d_Tresca_compute_error(FrictionContactProblem* problem, double* z, double* w,
                              double tolerance, SolverOptions* options, double* error);

#if defined(__cplusplus)
}
#endif

#endif
