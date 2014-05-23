/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#ifndef VariationalInequality_compute_error_H
#define VariationalInequality_compute_error_H

/*!\file VariationalInequality_computeError.h
  \brief functions related to error computation for friction-contact problems

  \author Vincent Acary, 26/05/2008

*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation for a VI problem. This function requires dWork to point to
   * at least 2*n double of allocated memory or it malloc this memory
      \param problem the structure which defines the VI problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param options
      \param[in,out] error value
      \return 0 if ok
   */
  int variationalInequality_computeError(VariationalInequality* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);

  /** Error computation for a box VI problem, that is\f$ \Pi_box(x-F(x)) - x\f$
      \param problem the structure which defines the VI problem
      \param x vector
      \param F vector
      \param tolerance value for error computation
      \param[in,out] error value
      \return 0 if ok
   */
  int variationalInequality_compute_error_box(
  VariationalInequality* problem,
  double *z , double *w, double tolerance, double* error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
