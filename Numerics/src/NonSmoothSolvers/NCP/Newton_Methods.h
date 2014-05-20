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

#ifndef Newton_Methods_H
#define Newton_Methods_H

#include "SolverOptions.h"

typedef void (*compute_F_ptr) (void* data_opaque, double* z, double* w);
typedef void (*compute_F_merit_ptr) (void* data_opaque, double* z, double* F, double* F_merit);

/** \struct _functions_FBLSA Newton_Methods.h
 * Struct holding the necessary pointer to functions
 */
typedef struct _functions_FBLSA {
  compute_F_ptr compute_F; /**< function to evaluate w = F(z) */
  compute_F_merit_ptr compute_F_merit; /**< function to evaluate F_merit(z) (e.g. F_FB, F_{min}, ...) */
  void (*compute_H)(void* data_opaque, double* z, double* w, double* workV1, double* workV2, double* H); /**< function to get an element H of T */
  void (*compute_error)(void* data_opaque, double* z, double* w, double* nabla_theta, double tol, double* err); /**< function to compute the error */
  void (*compute_RHS_desc)(void* data_opaque, double* z, double* w, double* F_desc); /**< function to evaluate F_desc(z) (e.g. F_FB, F_{min}, ...), optional */
  void (*compute_H_desc)(void* data_opaque, double* z, double* w, double* workV1, double* workV2, double* H_desc); /**< function to get an element H_desc of T_desc, optional */
} functions_FBLSA;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** NCP Solver based on FB with a line search
   * \param n size of the problem
   * \param z variable
   * \param w value of F(z)
   * \param info solver-specific values
   * \param data opaque problem definition
   * \param options options for this solver
   * \param functions struct of functions to compute F, H and the error
   */
  void newton_FBLSA(unsigned int n, double *z, double *w, int *info, void* data, SolverOptions* options, functions_FBLSA* functions);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
