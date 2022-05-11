/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef FRICTIONCONTACT3D_NONSMOOTH_NEWTON_SOLVERS_H
#define FRICTIONCONTACT3D_NONSMOOTH_NEWTON_SOLVERS_H

/*!\file fc3d_nonsmooth_Newton_solvers.h

  \brief Non Smooth Newton Solvers for FC3D
  [...]
 */

#include "NumericsFwd.h"  // for NumericsMatrix, FrictionContactProblem, Sol...

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




void computeAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB);
#endif
