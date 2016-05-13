/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

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
 * \author Olivier Huber
 */

#include "Newton_Methods.h"

/** \struct search_data ArmijoSearch.h
 * Struct to hold together the data needed by the search
 */
typedef struct {
  compute_F_ptr compute_F; /**< function to compute F(z) */
  compute_F_merit_ptr compute_F_merit; /**< function to compute F_merit(z) */
  double* z; /**< z vector */
  double* zc; /**< candidate z vector */
  double* F; /**< value of F(z) */
  double* F_merit; /**< value of F_merit(z) */
  double* desc_dir; /**< descent direction */
  double alpha0; /**< starting value for alpha */
  double alpha_min; /**< lower bound for alpha */
  void* data; /**< opaque pointer for extra data (needed for function call) */
  void* nm_ref_data; /**< data for the update rule */
  unsigned searchtype; /**< type of search: LINESEARCH or ARCSEARCH */
  void* set; /**< set on which the solution has to belong; only used for arcsearch */
  double sigma; /**< sigma value, used only for ARCSEARCH */
} search_data;


/** \struct nm_ref_struct ArmijoSearch.h
 * Struct used for the non-monotone search
 */
typedef struct {
  int type; /**< 0 if false, otherwise use a nonmonotone search. The integer value gives the update rule for the merit value ``threshold'' */
  int M; /**< maximum number of previous values of the merit function stored*/
  int m; /**< number of previous values of the merit function stored*/
  double* previous_thetas; /**< set of previous values of the merit function */
} nm_ref_struct;

enum { NM_LS_DISABLE, NM_LS_MAX, NM_LS_MEAN, NM_LS_ZHANG_HAGER };

enum SEARCH_TYPE { LINESEARCH, ARCSEARCH, BACKWARD_PATHSEARCH };

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Armijo linesearch; this version compute and update the reference value
   * and calls search_Armijo_standalone()
   * \param n size of the problem
   * \param theta current value of the merit function
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */
  double linesearch_Armijo2(int n, double theta, double preRHS, search_data* ls_data);

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
  double arcsearch_Armijo2(int n, double theta, double preRHS, search_data* ls_data);

  /** Armijo (non-monotone) search, standalone version: it does not compute
   * the reference value, it is expected as argument (theta)
   * \param n size of the problem
   * \param theta reference value for the acceptance test
   * \param preRHS pre-computed value for the acceptance test
   * \param ls_data necessary data for the search algorithm
   * \return the coefficient alpha
   */
  double search_Armijo_standalone(int n, double* theta, double preRHS, search_data* ls_data);

  /** Set the update type for the non-monotone line search
   * \param nm_ref_data the struct containing the data relative to the
   * reference value used in the line search
   * \param type for the update
   */
  static inline void set_nonmonotone_type(void* nm_ref_data, int type)
  {
    ((nm_ref_struct*) nm_ref_data)->type = type;
  }

  /** Get the update type for the non-monotone line search
   * \param nm_ref_data the struct containing the data relative to the
   * reference value used in the line search
   * \return the type for the update
   */
  static inline int get_nonmonotone_type(void* nm_ref_data)
  {
    return ((nm_ref_struct*) nm_ref_data)->type;
  }

  /** update the reference value for the non-monotone line search
   * \param nm_ref_data the struct containing the data for the update
   * \param cur_merit new value of the merit function for the update
   */
  void update_non_monotone_ref(void* nm_ref_data, double cur_merit);

  /** compute the reference value
   * \param nm_ref_data the struct containing the data relative to the
   * reference value used in the line search
   * \param[in,out] theta_ref on input the current value of the merit function;
   * on output, the new reference value
   */
  void get_non_monotone_ref(void* nm_ref_data, double* theta_ref);

  /** fill the data struct for non-monotone search
   * \param nm_ref_data the structure to fill
   * \param iparam the set of parameter from the SolverOption struct
   */
  void fill_nm_data(nm_ref_struct* nm_ref_data, int* iparam);

  /** free the allocated memory for the non-monotone search
   * \param nm_ref_data the structure holding the allocated memory
   */
  void free_nm_data(nm_ref_struct* nm_ref_data);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
