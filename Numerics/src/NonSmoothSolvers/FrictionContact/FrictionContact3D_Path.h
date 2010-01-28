/* Siconos-Numerics, Copyright INRIA 2005-2010.
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

#ifdef __cplusplus
extern "C" {
#endif
  /** Initialize friction-contact 3D Path solver
      \param dim. of the global problem
      \param matrix M of the global problem
      \param vector q of the global problem
      \param vector of the friction coefficients
      \param  SolverOptions * options of the solver
   */
  void frictionContact3D_Path_initialize(int, const NumericsMatrix*const, const double*const, const double*const, SolverOptions *);

  /** solve friction-contact 3D problem with Path
      \param number (position in global matrix) of the considered contact
      \param dim. of the global problem
      \param global reaction (only the block corresponding to the current contact will be modified,
      \param vector of int parameters (max iteration numnber ...)
      \param  SolverOptions * options of the solver
   */
  void frictionContact3D_Path_solve(int, int, double*, SolverOptions*);

  /** free memory for friction contact 3D Path solver */
  void frictionContact3D_Path_free();

  /** solve friction-contact 3D problem with Path
      \param dimension of the global problem
      \param velocity vector (in-out)
      \param global reaction vector
      \param output error
  */
  void frictionContact3D_Path_computeError(int, double*, double*, double *);

#ifdef __cplusplus
}
#endif

#endif
