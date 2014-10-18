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
#include "AVI_Solvers.h"

/*!\file avi_caoferris.h 
 *  \brief Subroutines for the solver by Cao and Ferris
 */

#ifndef AVI_CAOFERRIS_H
#define AVI_CAOFERRIS_H

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** stage 3 of the Cao-Ferris algorithm
   * \param problem struct formalizing the AVI
   * \param u vector for the basic variables
   * \param s vector the non-basic variables
   * \param size_x dimension of the solution variable
   * \param A set of active constraints
   * \param options struct used to define the solver(s) and its (their)
   * parameters
   * \return 0 if success, 1 if failure
   * \author Olivier Huber
   */
  int avi_caoferris_stage3(AffineVariationalInequalities* problem, double* u , double* s, unsigned int size_x, unsigned int* A, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
