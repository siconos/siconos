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

#ifndef NCP_NEWTON_FB
#define NCP_NEWTON_FB

/*!\file NCP_newton_FBLSA.h
 * \brief functions for solving NCP using a Newton method on the nonsmooth
 * equation reformulation*/

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Compute F(z) = Mz + q
   * \param data_opaque a LinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[out] w basic variable (result)
   */
  void FB_compute_F_ncp(void* data_opaque, double* z, double* w);

  /** Compute an element of JacF_FB, see Facchinei--Pang p. 808
   * \param data_opaque a NonlinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[in] w basic variable
   * \param workV1 work vector which contains "z"
   * \param workV2 work vector
   * \param[out] H an element of JacF_FB
   */
  void FB_compute_H_ncp(void* data_opaque, double* z, double* w, double* workV1, double* workV2, NumericsMatrix* H);

  /** Compute the error for termination, here ncp_compute_error
   * \param data_opaque a NonlinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[in] w basic variable
   * \param notused not used here
   * \param[in] tol the tolerance
   * \param[out] err the error on the LCP (not FB)
   */
  void FB_compute_error_ncp(void* data_opaque, double* z, double* w, double* notused, double tol, double* err);

  /** Compute F_FB : \f${F_FB}_i = \sqrt(z_i^2 + F_i^2) - (z_i + F_i)\f$
   * \param data_opaque a NonlinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[in] F basic variable
   * \param[out] F_FB value of the function
   */
  void ncp_FB(void* data_opaque, double* z, double* F, double* F_FB);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
