/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#ifndef grfc3d_compute_error_H
#define grfc3d_compute_error_H

/*!\file grfc3d_compute_error.h
  \brief functions related to error computation for rolling friction-contact problems

*/

#include "NumericsFwd.h"    // for GlobalRollingFrictionContactProblem, SolverOptions
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation (using the normal map residual) for rolling friction-contact 3D problem
      \param problem the structure which defines the rolling friction-contact problem
      \param u vector
      \param r vector
      \param tolerance value for error computation
      \param options
      \param norm of a vector (problem->q) for relative error
      \param[in,out] error value
      \return 0 if ok
   */
  int grfc3d_compute_error(GlobalRollingFrictionContactProblem* problem,
                          double *reaction , double *velocity,
                          double* globalVelocity, double tolerance,
                          double * error, int problemIsNotConvex);

  /** Error computation (using the normal map residual) for one rolling friction-contact 3D problem
      \param r the reaction force
      \param u the local velocity
      \param mu coeficient of friction
      \param mur coeficient of rolling
      \param worktmp work vector
      \param[in,out] error value
   */
  void grfc3d_unitary_compute_and_add_error(double* r, double* u, double mu,
                                            double mur, double * error, double * worktmp,
                                            int problemIsNotConvex);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
