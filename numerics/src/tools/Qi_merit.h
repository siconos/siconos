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

#ifndef QI_MERIT_H
#define QI_MERIT_H

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#ifdef __cplusplus
#define restrict __restrict
#endif

/*!\file Qi_merit.h
  \brief functions related to the Qi C-functions used as a merit function for box VI problems

  Reference Facchinei--Pang pp.869 - 877.

*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Evaluates the C function for a box-constrained VI
   * \param n size of the problem
   * \param[in] x box-constrained variable of the VI
   * \param[out] F value of the function
   * \param[out] Fbox value of the function
   * \param[in] lb lower bounds, that is lb <= x
   * \param[in] ub upper bounds, that is ub >= x
   * */
  void phi_Qi(int n, double* restrict x, double* restrict F, double* restrict Fbox, double* restrict lb, double* restrict ub);

  /** Evaluates the Jacobian of the C function for a box-constrained VI
   * \param n size of the problem
   * \param[in] x box-constrained variable of the VI
   * \param[out] Fbox value of the function
   * \param workV1 work vector
   * \param workV2 work vector
   * \param[in] nabla_F gradient of the C-function
   * \param[in] lb lower bounds, that is lb <= x
   * \param[in] ub upper bounds, that is ub >= x
   * \param[out] H an element of the Jacobian
   * */
  void Jac_F_Qi(int n, double* restrict x, double* restrict Fbox, double* restrict workV1, double* restrict workV2, NumericsMatrix* restrict nabla_F, double* restrict lb, double* restrict ub, NumericsMatrix* restrict H);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
