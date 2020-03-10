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

/** \file pivot-utils.h 
 * \brief some functions used in one or more solver, mainly related to Lemke's
 * algorithm
 */

#ifndef PIVOT_UTILS_H
#define PIVOT_UTILS_H

#include "NumericsFwd.h"    // for NumericsMatrix
#include "lumod_wrapper.h"  // for SN_lumod_dense_data
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#ifdef __cplusplus
#undef restrict
#include <sys/cdefs.h>      // for __restrict
#define restrict __restrict
#endif

#define PIVOT_PATHSEARCH_SUCCESS -2

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

/** find the first leaving variable for Lemke's algorithm
 * \param mat the matrix
 * \param dim the dimension of the problem
 * \return the index of the blocking variable
 */
  int pivot_init_lemke(double* mat, unsigned int dim);

/** find the first leaving variable for the pathsearch algorithm
 * \param dim the dimension of the problem
 * \param mat the matrix
 * \param[out] t_indx the index of t in the basis (needed since block < 0 is
 * possible and meaningful)
 * \return the index of the blocking (leaving) variable
 */
  int pivot_init_pathsearch(unsigned dim, double* mat, unsigned* t_indx);

/** find the leaving variable for Lemke's algorithm
 * \param mat the matrix; we only use the columns of basic variable and the
 * covering vector of the entering variable
 * \param dim dimension of the problem
 * \param drive the driving variable
 * \param aux_index index of auxillary variable in the current basis
 * \return the leaving (or blocking) variable
 */
  int pivot_selection_lemke(double* mat, unsigned dim, unsigned drive, unsigned aux_index);

  /** find the leaving variable for Lemke's algorithm
 * \param n dimension of the problem
 * pivots)
 * \param col_drive column of the driving variable (contains all the possible
 * pivots)
 * \param q_tilde solution of the linear system
 * \param lexico_mat matrix for the lexico ordering
 * \param aux_indx index of auxillary variable in the current basis
 * \param lexico_tol the tolerance on the lexicographic comparison
 * \return the leaving (or blocking) variable
 */
  int pivot_selection_lemke2(unsigned n, double* col_drive, double* q_tilde, double* lexico_mat, unsigned aux_indx, double lexico_tol);

  /** find the leaving variable for Lemke's algorithm
 * \param n dimension of the problem
 * \param col_drive column of the driving variable (contains all the possible
 * pivots)
 * \param q_tilde solution of the linear system
 * \param lexico_col column for the lexico ordering
 * \param basis current basis
 * \param candidate_indx array for storing the possible candidate indexes
 * \param lumod_data data for the BLU update
 * \param aux_indx index of auxillary variable in the current basis
 * \param lexico_tol the tolerance on the lexicographic comparison
 * \return the leaving (or blocking) variable
 */
int pivot_selection_lemke3(unsigned n, double* col_drive, double* q_tilde, double* lexico_col, unsigned* basis, unsigned* candidate_indx, SN_lumod_dense_data* lumod_data, unsigned aux_indx, double lexico_tol);

/** find the leaving variable in a path search procedure. The code is almost
 * the same as in pivot_selection_lemke, except that we also check if it is
 * possible to set t = 1.
 * \param mat the matrix (we only use the columns of basic variable and the
 * covering vector of the entering variable
 * \param dim dimension of the problem
 * \param drive the driving or entering variable
 * \param t_indx the index of the t variable in the basis
 * \return the leaving (or blocking) variable
 */
  int pivot_selection_pathsearch(double* mat, unsigned dim, unsigned drive, unsigned t_indx);

/** Initialize the matrix for Lemke's algorithm as [ q | Id | -d | -M ] with
 * d_i = 1 if i < size_x if the argument d is NULL. Otherwise take d from the
 * arguments
 * \param mat column-major matrix to fill
 * \param M M matrix from the LCP formulation
 * \param dim dimension of the problem
 * \param size_x for the classical Lemke's algorithm, this is set to dim. Its
 * value is different when we solve an AVI using Cao and Ferris method
 * \param q vector from the LCP formulation
 * \param d covering vector for the auxiliary variable, possibly NULL
 */
void init_M_lemke(double* mat, double* M, unsigned int dim, unsigned int size_x, double* q, double* d);

/** Do the pivot <block, drive>, driftless version
 * \param mat the matrix to update
 * \param dim the number of rows
 * \param dim2 the number of columns
 * \param block the blocking or leaving variable
 * \param drive the driving or entering variable
 */
void do_pivot_driftless(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive);

/** Do the pivot <block, drive>, driftless version
 * \param mat the matrix to update
 * \param dim the number of rows
 * \param dim2 the number of columns
 * \param block the blocking or leaving variable
 * \param drive the driving or entering variable
 */
void do_pivot_driftless2(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive);

/** Do the pivot <block, drive>, standart version
 * \param mat the matrix to update
 * \param dim the number of rows
 * \param dim2 the number of columns
 * \param block the blocking or leaving variable
 * \param drive the driving or entering variable
 */
void do_pivot(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive);

/** Do the pivot <block, drive> with block-LU updates
 * \param lumod_data lumod data
 * \param M the LCP matrix
 * \param q_tilde the modified q vector: it is the current of the variables currently in the basis
 * \param lexico_mat matrix for the lexicographic ordering
 * \param col_drive column of the driving variable
 * \param col_tilde Solution to H x = col
 * \param basis current basis
 * \param block the blocking or leaving variable
 * \param drive the driving or entering variable
 */
void do_pivot_lumod(SN_lumod_dense_data* lumod_data, NumericsMatrix* M, double* q_tilde, double* lexico_mat, double* col_drive, double* col_tilde, unsigned* basis, unsigned block, unsigned drive);

  /** print a diagnostic of the info value
   * \param info the value given by the algorithm
   */
  void lcp_pivot_diagnose_info(int info);

  int pivot_selection_bard(double* mat, unsigned int dim);
  int pivot_selection_least_index(double* mat, unsigned int dim);
  void init_M_bard(double* restrict mat, double* restrict M, unsigned int dim, double* restrict q);
  void init_M_least_index(double* restrict mat, double* restrict M, unsigned int dim, double* restrict q);
  int init_M_lemke_warm_start(int n, double* restrict u, double* restrict mat, double* restrict M, double* restrict q, int* restrict basis, double* restrict cov_vec);

  const char* basis_to_name(unsigned nb, unsigned n);
  unsigned basis_to_number(unsigned nb, unsigned n);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
