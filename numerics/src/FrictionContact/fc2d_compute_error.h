/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "SiconosConfig.h"
#include "FrictionContactProblem.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation for friction-contact 3D problem
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param[in,out] error value
      \return 0 if ok
   */
  int fc2d_compute_error(FrictionContactProblem* problem, double *z , double *w, double tolerance, double norm, double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
