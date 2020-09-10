/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#ifndef rolling_fc3d_compute_error_H
#define rolling_fc3d_compute_error_H

/*!\file rolling_fc3d_compute_error.h
  \brief functions related to error computation for friction-contact problems

*/

#include "NumericsFwd.h"  // for RollingFrictionContactProblem, SolverOptions
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation (using the normal map residual) for rolling friction-contact 2D problem
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param options
      \param norm norm of a vector (problem->q) for relative error
      \param[in,out] error value
      \return 0 if ok
   */
  int rolling_fc2d_compute_error(RollingFrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double norm, double * error);

  /** Error computation (using the normal map residual) for one rolling friction-contact 2D problem
      \param r the reaction force
      \param u the local velocity
      \param mu coeficient of friction
      \param worktmp work vector
      \param[in,out] error value
   */
  void rolling_fc2d_unitary_compute_and_add_error(double r[3] , double u[3], double mu, double mur, double * error, double * worktmp);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
