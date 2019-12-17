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

#ifndef VariationalInequality_compute_error_H
#define VariationalInequality_compute_error_H

/*!\file VariationalInequality_computeError.h
  \brief functions related to error computation for friction-contact problems

*/

#include "NumericsFwd.h"  // for VariationalInequality, SolverOptions
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation for a VI problem. This function requires dWork to point to
   * at least 2*n double of allocated memory or it malloc this memory
      \param problem the structure which defines the VI problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param options solver options
      \param[in,out] error value
      \return 0 if ok
   */
  int variationalInequality_computeError(VariationalInequality* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);

  /** Error computation for a box VI problem, that is\f$ \Pi_box(x-F(x)) - x\f$
      \param problem the structure which defines the VI problem
      \param x vector
      \param F vector
      \param tolerance value for error computation
      \param[in,out] error value
      \return 0 if ok
   */
  int variationalInequality_compute_error_box(VariationalInequality* problem, double* restrict x , double* restrict F, double tolerance, double* restrict error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
