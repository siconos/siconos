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
#ifndef FRICTIONCONTACT3D_local_nonsmooth_Newton_solvers_H
#define FRICTIONCONTACT3D_local_nonsmooth_Newton_solvers_H

/*!\file fc3d_local_nonsmooth_Newton_solvers.h
  \brief Typedef and functions declarations related to Newton solver for 3 dimension frictional contact problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for a specific contact
  - solve: solve the local problem
  - free

*/
#include "NumericsMatrix.h"
#include "SolverOptions.h"
#include "FrictionContactProblem.h"
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

typedef void (*computeNonsmoothFunction)(double *, double * , double , double * , double *, double *, double *);

  /** initialize friction-contact 3D Newton solver
   * \param problem to solve
   * \param localproblem to solve
   * \param options of the solver
   */
  void fc3d_local_nonsmooth_Newton_solvers_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * options);

  /** solve friction-contact 3D problem with Newton
   * \param localproblem to solve
   * \param options of the solver
   * \return 0 iff successful.
   */
  int fc3d_local_nonsmooth_Newton_solvers_solve(FrictionContactProblem* localproblem, double*, SolverOptions * options);

  /** free memory for friction contact 3D Newton solver
   * \param localproblem for freeing matrix0
   */
  void fc3d_local_nonsmooth_Newton_solvers_free(FrictionContactProblem* localproblem);

  /** compute error for friction-contact 3D problem with Newton
   *  \param dimension of the global problem
   *   \param[in,out] velocity vector (\warning in-out parameter )
   *   \param reaction global reaction vector
   *   \param output_error
   */
  void fc3d_local_nonsmooth_Newton_solvers_computeError(int dimension, double* velocity, double*reaction, double * output_error);

  /** Update friction-contact 3D problem: formalize local problem for one contact
      \param problem the global problem to solve
      \param localproblem the local problem to solve
      \param number (position in global matrix) of the considered contact
      \param reaction global reaction (only the block corresponding to the
      current contact will be modified
      \param options of the solver

      the rest is used to formalize the local problem)
  */
  void fc3d_local_nonsmooth_Newton_AC_update(int number, FrictionContactProblem* problem, FrictionContactProblem* localproblem ,
                                   double * reaction, SolverOptions* options);

  int fc3d_local_nonsmooth_Newton_solvers_solve_direct(FrictionContactProblem* localproblem,
                                 double * R, int *iparam, double *dparam);

  int fc3d_local_nonsmooth_Newton_solvers_solve_damped(FrictionContactProblem* localproblem,
                                       double * R, int *iparam, double *dparam);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
