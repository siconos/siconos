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

#ifndef ARMIJOSEARCH_H
#define ARMIJOSEARCH_H

/*!\file ArmijoSearch.h
 * \brief Armijo type search (linesearch and arcsearch), possibly non-monotone
 *
 * Implementation of the Armijo type search. The linesearch version is the most
 * classical one. The arcsearch assumes that the direction of descent is the
 * gradient of the merit function.
 *
 */

#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep
#include "line_search.h"

/** \struct armijo_extra_params ArmijoSearch.h
 * Struct to hold together the extra parameters needed by the Armijo line search
 */
typedef struct {
  double gamma; /**< Value of the slope coefficient*/
} armijo_extra_params;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Armijo (non-monotone) search, standalone version: it does not compute
   * the reference value, it is expected as argument (theta)
   * \param n size of the problem
   * \param theta reference value for the acceptance test
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */
  double search_Armijo_standalone(int n, double* theta, double preRHS, search_data* ls_data);

  /** Armijo linesearch; this version compute and update the reference value
   * and calls search_Armijo_standalone()
   * \param n size of the problem
   * \param theta current value of the merit function
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */
  static inline double linesearch_Armijo2(int n, double theta, double preRHS, search_data* ls_data)
  {
    return line_search_generic(n, theta, preRHS, ls_data, LINESEARCH, &search_Armijo_standalone);
  }

  /** Armijo arcsearch; this version compute and update the reference value
   * and calls search_Armijo_standalone().
   * \warning this function can be used only if the descent direction is the
   * gradient of the merit function.
   * \param n size of the problem
   * \param theta current value of the merit function
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */

  static inline double arcsearch_Armijo2(int n, double theta, double preRHS, search_data* ls_data)
  {
    return line_search_generic(n, theta, preRHS, ls_data, ARCSEARCH, &search_Armijo_standalone);
  }

  /** Initialize parameters to a default value
   * \param p parameters to set
   */
  void search_Armijo_params_init(armijo_extra_params* p);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
