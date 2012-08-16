/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#ifndef FRICTIONCONTACT3DPath_H
#define FRICTIONCONTACT3DPath_H

/*!\file FrictionContact3D_Path.h
  \brief Typedef and functions declarations related to NCP-Path solver for 3 dimension frictional contact problems.
  \author Franck Perignon

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for a specific contact
  - solve: solve the local problem
  - free

*/
#include "SparseBlockMatrix.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Initialize friction-contact 3D Path solver
   * \param problem to solve
   * \param localproblem to solve
   * \param localsolver_options of the solver
   */
  void frictionContact3D_Path_initialize(FrictionContactProblem * problem , FrictionContactProblem * localproblem, SolverOptions * localsolver_options);

  /** solve friction-contact 3D problem with Path
   * \param localproblem to solve
   * \param reaction
   * \param options of the solver
   */
  int frictionContact3D_Path_solve(FrictionContactProblem * localproblem , double* reaction, SolverOptions* options);

  /** free memory for friction contact 3D Path solver */
  void frictionContact3D_Path_free();

  /**  compute error for  friction-contact 3D problem with Path
   * \param dimension of the global problem
   * \param[in,out] velocity vector (\warning in-out parameter )
   * \param reaction vector
   * \param output_error
   */
  void frictionContact3D_Path_computeError(int dimension, double* velocity, double* reaction, double * output_error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
