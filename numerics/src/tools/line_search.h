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

#ifndef LINESEARCH_H
#define LINESEARCH_H

/*!\file line_search.h
 * \brief Basic structures for line-search (and arcsearch) methods
 *
 */

#include "Newton_methods.h"
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

/** \struct search_data line_search.h
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
  double alpha_min; /**< minimum value of alpha*/
  void* data; /**< opaque pointer for extra data (needed for function call) */
  void* nm_ref_data; /**< data for the update rule */
  unsigned searchtype; /**< type of search: LINESEARCH or ARCSEARCH */
  void* set; /**< set on which the solution has to belong; only used for arcsearch */
  double sigma; /**< sigma value, used only for ARCSEARCH */
  void* extra_params; /**< extra parameters for some line search algorithm */
} search_data;


/** \struct nm_ref_struct line_search.h
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

enum LSA_ALGO { SICONOS_LSA_ARMIJO, SICONOS_LSA_GOLDSTEIN };

typedef double (*sn_ls_fn)(int n, double* theta, double preRHS, search_data* ls_data);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

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

  /** Generic call for a linesearch (or arcsearch). Handles the update of the
   * update of the non-monotone data
   * \param n size of the variable
   * \param theta current value of the merit function
   * \param preRHS value used for the comparison
   * \param ls_data line search data
   * \param searchtype type of search: linesearch or arcsearch
   * \param ls_fn function to call
   * \return the value of tau, NAN if the search failed
   */
  double line_search_generic(int n, double theta, double preRHS, search_data* ls_data, unsigned searchtype, sn_ls_fn ls_fn);

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

  /** Reset the storage for the non-monotone search
   * \param nm_ref_data the structure
   */
  void zero_nm_data(nm_ref_struct* nm_ref_data);

  /** free the allocated memory for the non-monotone search
   * \param nm_ref_data the structure holding the allocated memory
   */
  void free_nm_data(nm_ref_struct* nm_ref_data);

  /** free the allocated memory for the linesearch method
   * \param ls_data the struct
   */
  void free_ls_data(search_data* ls_data);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

