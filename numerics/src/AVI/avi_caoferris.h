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

/*!\file avi_caoferris.h 
 *  \brief Subroutines for the solver by Cao and Ferris
 */

#ifndef AVI_CAOFERRIS_H
#define AVI_CAOFERRIS_H

#include "AVI_Solvers.h"
#include "LinearComplementarityProblem.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** stage 3 of the Cao-Ferris algorithm
   * \param problem struct formalizing the AVI
   * \param u vector for the basic variables
   * \param s vector the non-basic variables
   * \param d the covering vector for Lemke
   * \param size_x dimension of the solution variable
   * \param A set of active constraints
   * \param options struct used to define the solver(s) and its (their) parameters
   * \return 0 if success, 1 if failure
   */
  int avi_caoferris_stage3(LinearComplementarityProblem* problem, double* u, double* s, double* d, unsigned size_x, unsigned* A, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
