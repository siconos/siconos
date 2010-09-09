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
#ifndef FRICTIONCONTACT3DUNITARY_ENUMERATIVE_H
#define FRICTIONCONTACT3DUNITARY_ENUMERATIVE_H

/*!\file FrictionContact3D_Path.h
  \brief Typedef and functions declarations related to the quartic solver for 3 dimension frictional contact problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for a specific contact
  - solve: solve the local problem
  - free

*/
#include "SparseBlockMatrix.h"
#include "SolverOptions.h"

#ifdef __cplusplus
extern "C"
{
#endif
  void frictionContact3D_unitary_enumerative_free(FrictionContactProblem* problem);
  void frictionContact3D_unitary_enumerative_initialize(FrictionContactProblem* problem);
  int frictionContact3D_unitary_enumerative_solve(FrictionContactProblem* problem, double * reaction, SolverOptions* options);
  int frictionContact3D_unitary_enumerative(FrictionContactProblem* problem, double * reaction, double * velocity, int *info, SolverOptions* options);
  int frictionContact3D_unitary_enumerative_setDefaultSolverOptions(SolverOptions* options);
#ifdef __cplusplus
}
#endif

#endif
