/* Siconos-Numerics, Copyright INRIA 2005-2014.
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

#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "Newton_Methods.h"

/** \struct _linesearch_data LineSearch.h
 * Struct to hold together the data needed by the linesearch
 */
typedef struct _linesearch_data {
  compute_F_ptr compute_F; /** function to compute F(z) */
  compute_F_merit_ptr compute_F_merit; /** function to compute F_merit(z) */
  double* z; /**< z vector */
  double* zc; /**< candidate z vector */
  double* F; /**< value of F(z) */
  double* F_merit; /**< value of F_merit(z) */
  double* desc_dir; /**< descent direction */
  void* data; /**< opaque pointer for extra data (needed for function call*/
  int nonmonotone; /**< 0 if false, otherwise use os nonmonotone linesearch. The integer value gives the update rule for the merit value ``threshold'' */
  int M; /**< maximum number of previous values of the merit function stored*/
  int m; /**< number of previous values of the merit function stored*/
  double* previous_thetas; /**< set of previous values of the merit function */
} linesearch_data;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Armijo linesearch
   * \param n size of the problem
   * \param theta current value of the merit function
   * \param ls_data necessary data for the linesearch algorithm
   */
  double linesearch_Armijo2(int n, double theta, double preRHS, linesearch_data* ls_data);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
