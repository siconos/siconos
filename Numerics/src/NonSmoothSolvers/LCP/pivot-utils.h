/* Siconos-Numerics, Copyright INRIA 2005-2014
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/** \file pivot-utils.h 
 * \brief some functions used in one or more solver, mainly related to Lemke's
 * algorithm
 */

#ifndef PIVOT_UTILS_H
#define PIVOT_UTILS_H

#include "NumericsConfig.h"

#ifdef __cplusplus
#undef restrict
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
 * \param drive the driving or entering variable
 * \return the leaving (or blocking) variable
 */
  int pivot_selection_lemke(double* mat, unsigned int dim, unsigned int drive);

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
  int pivot_selection_pathsearch(double* mat, unsigned int dim, unsigned int drive, unsigned int t_indx);

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
void init_M_lemke(double* restrict mat, double* restrict M, unsigned int dim, unsigned int size_x, double* restrict q, double* restrict d);

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

  /** print a diagnostic of the info value
   * \param info the value given by the algorithm
   */
  void lcp_pivot_diagnose_info(int info);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
