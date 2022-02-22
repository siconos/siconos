/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifndef GOLDSTEINSEARCH_H
#define GOLDSTEINSEARCH_H

/*!\file GoldsteinSearch.h
 * \brief Goldstein type search (linesearch and arcsearch), possibly non-monotone
 *
 * Implementation of the Goldstein type search. The linesearch version is the most
 * classical one. The arcsearch assumes that the direction of descent is the
 * gradient of the merit function.
 *
 */

#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep
#include <stddef.h>       // for size_t
#include "line_search.h"  // for line_search_generic, search_data, ARCSEARCH

/** \struct goldstein_extra_params GoldsteinSearch.h
 * Struct to hold together the extra parameters needed by the Goldstein line search
 */
typedef struct {
  size_t iter_max; /**< maximum number of iterations */
  double c; /**< Value of the slope coefficient*/
  double alpha_max; /**< maximum value of alpha*/
} goldstein_extra_params;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Goldstein (non-monotone) search, standalone version: it does not compute
   * the reference value, it is expected as argument (theta)
   * \param n size of the problem
   * \param theta reference value for the acceptance test
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */
  double search_Goldstein_standalone(int n, double* theta, double preRHS, search_data* ls_data);

  /** Goldstein linesearch; this version compute and update the reference value
   * and calls search_Goldstein_standalone()
   * \param n size of the problem
   * \param theta current value of the merit function
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */
  static inline double linesearch_Goldstein2(int n, double theta, double preRHS, search_data* ls_data)
  {
    return line_search_generic(n, theta, preRHS, ls_data, LINESEARCH, &search_Goldstein_standalone);
  }

  /** Goldstein arcsearch; this version compute and update the reference value
   * and calls search_Goldstein_standalone().
   * \warning this function can be used only if the descent direction is the
   * gradient of the merit function.
   * \param n size of the problem
   * \param theta current value of the merit function
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */
  static inline double arcsearch_Goldstein2(int n, double theta, double preRHS, search_data* ls_data)
  {
    return line_search_generic(n, theta, preRHS, ls_data, ARCSEARCH, &search_Goldstein_standalone);
  }

  /** Initialize parameters to a default value
   * \param p parameters to set
   */
  void search_Goldstein_params_init(goldstein_extra_params* p);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
