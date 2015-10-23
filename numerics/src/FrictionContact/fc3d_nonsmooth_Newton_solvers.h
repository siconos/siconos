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
#ifndef FRICTIONCONTACT3D_NONSMOOTH_NEWTON_SOLVERS_H
#define FRICTIONCONTACT3D_NONSMOOTH_NEWTON_SOLVERS_H

/*!\file fc3d_nonsmooth_Newton_solvers.h

  \brief Non Smooth Newton Solvers for FC3D
  [...]
 */

#include "FrictionContactProblem.h"
#include "SolverOptions.h"

/* Notes:
 - check Olivier Newton_method (only dense)
 x only sparse storage. (fixed now)
 - void parameter *data is unused for the moment.
 - problem_size, mu and rho should be moved outside parameters (inside *data)
*/


/** The nonsmooth function signature.
[...]
 */
typedef void (*fc3d_nonsmooth_Newton_solversFunPtr)
(void* data,
 unsigned int problem_size,
 double* reaction,
 double* velocity,
 double* mu,
 double* rho,
 double* F,
 double* A,
 double* B);

/** The nonsmooth function signature for a 3x3 block.
 */
typedef void (*FrictionContactNSFun3x3Ptr)(double* reaction,
                                           double* velocity,
                                           double mu,
                                           double* rho,
                                           double* F,
                                           double* A,
                                           double* B);

/** The nonsmooth equation structure.
[...]
 */
typedef struct
{
  FrictionContactProblem* problem;
  void* data;
  fc3d_nonsmooth_Newton_solversFunPtr function;
} fc3d_nonsmooth_Newton_solvers;

/** Solve the equation. The only implemented method is a nonsmooth
    Newton method with a Goldstein Price or a FBLSA line search.
    Linear solver choice and line search are specified in
    SolverOptions parameter.
    \param equation the nonsmooth equation.
    \param reaction the reaction guess as input and the solution as output.
    \param velocity the velocity guess as input and the solution as output.
    \param info the return info. 0 success, 1 failure.
    \param options the SolverOptions parameter.
 */
void fc3d_nonsmooth_Newton_solvers_solve(fc3d_nonsmooth_Newton_solvers* equation,
                                      double* reaction,
                                      double* velocity,
                                      int* info,
                                      SolverOptions* options);




#endif
