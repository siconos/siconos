/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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



#ifndef SolverOptions_helpers_H
#define SolverOptions_helpers_H

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** set the tolerance for the solver
   * \param dparam the set of double parameters
   * \param tol the new tolerance
   */
  inline static void SN_set_tolerance(double* dparam, double tol) { dparam[0] = tol; };

  /** get the tolerance for the solver
   * \param dparam the set of double parameters
   * \return tol the tolerance
   */
  inline static double SN_get_tolerance(double* dparam) { return dparam[0]; };

  /** set the residual from the solver
   * \param dparam the set of double parameters
   * \param res the new tolerance
   */
  inline static void SN_set_residual(double* dparam, double res) { dparam[1] = res; };

  /** get the residual
   * \param dparam the set of double parameters
   * \return the residual
   */
  inline static double SN_get_residual(double* dparam) { return dparam[1]; };

  /** set the maximum number of iterations
   * \param iparam the set of integer parameters
   * \param max_iters the new tolerance
   */
  inline static void SN_set_max_iters(int* iparam, int max_iters) { iparam[0] = max_iters; };

  /** get the number of iterations 
   * \param iparam the set of integer parameters
   * \return the maximum number of iterations
   */
  inline static int SN_get_max_iters(int* iparam) { return iparam[1]; };

  /** set the number of iterations done in the solver
   * \param iparam the set of integer parameters
   * \param nb_iters the new tolerance
   */
  inline static void SN_set_nb_iters(int* iparam, int nb_iters) { iparam[1] = nb_iters; };

  /** get the number of iterations 
   * \param iparam the set of double parameters
   * \return the number of iterations done in the solver
   */
  inline static int SN_get_nb_iters(int* iparam) { return iparam[1]; };

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif


#endif
