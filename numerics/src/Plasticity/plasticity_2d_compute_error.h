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

#ifndef plasticity_2d_compute_error_H
#define plasticity_2d_compute_error_H

/*!\file plasticity_2d_compute_error.h
  \brief functions related to error computation for Mohr Coulomb 2D plasticity

*/

#include "NumericsFwd.h"  // for PlasticityProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/**
    Error computation (using the normal map residual) for Mohr Coulomb 2D plasticity

    \param problem the structure which defines the  Mohr Coulomb 2D plasticity
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param norm norm of a vector (problem->q) for relative error
    \param[in,out] error value
    \return 0 if ok
 */
int plasticity_2d_compute_error(PlasticityProblem *problem, double *z, double *w,
                                double tolerance, SolverOptions *options, double norm,
                                double *error);

/**
    Error computation (using the normal map residual) for one Mohr Coulomb 2D plasticity

    \param r the reaction force
    \param u the local velocity
    \param eta coefficient of the cone
    \param theta coefficient of the dilatency
    \param worktmp work vector
    \param[in,out] error value
 */
void plasticity_2d_unitary_compute_and_add_error(double *restrict r, double *restrict u,
                                                 double eta, double theta, double *error,
                                                 double *worktmp);

#if defined(__cplusplus)
}
#endif

#endif
