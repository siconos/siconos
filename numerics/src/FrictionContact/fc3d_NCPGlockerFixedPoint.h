/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef FRICTIONCONTACT3DNCPGlockerFixedPoint_H
#define FRICTIONCONTACT3DNCPGlockerFixedPoint_H

/*!\file fc3d_NCPGlockerFixedPoint.h
  \brief Typedef and functions declarations related to NCP-Fixed Point solver for 3 dimension frictional contact problems.
  \author Houari Khenous
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

  void F_GlockerFixedP(int sizeF, double* reaction, double* FVector, int up2Date);


  /** Initialize friction-contact 3D Fixed Point solver
   * \param problem to solve
   * \param localproblem to solve
   * \param localsolver_options of the solver
   */
  void fc3d_FixedP_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions * localsolver_options);

  /** solve friction-contact 3D problem with Fixed Point
      \param localproblem to solve
      \param reaction (only the block corresponding to the current contact will be modified,
      \param  options of the solver
      \return 0 iff successful
   */
  int fc3d_FixedP_solve(FrictionContactProblem * localproblem , double* reaction , SolverOptions * options);

  /** free memory for friction contact 3D Fixed Point solver */
  void fc3d_FixedP_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions * localsolver_option);

  /** compute error for friction-contact 3D problem with Fixed Point
   * \param dimension of the global problem
   * \param[in,out] velocity vector (\warning in-out parameter )
   * \param reaction vector
   * \param output_error
  */
  void fc3d_Path_computeError(int dimension, double* velocity, double* reaction, double * output_error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
