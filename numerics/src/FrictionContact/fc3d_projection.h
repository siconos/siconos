/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#ifndef FRICTIONCONTACT3DProjection_H
#define FRICTIONCONTACT3DProjection_H

/*!\file fc3d_projection.h
  \brief Typedef and functions declarations related to projection solver for 3 dimension frictional contact problems

  Each solver must have 4 functions in its interface:
  - initialize: link global static variables to the considered problem (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for a specific contact
  - solve: solve the local problem
  - free
  We consider a "global" (ie for several contacts) problem, used to initialize the static global variables.
  Then a "local" (ie for one contact => size = 3) problem is built (update function) and solved (solve function).

  Two different storages are available for M: dense and sparse block.

  \author INRIA Siconos Team

*/
#include "NumericsMatrix.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Initialize friction-contact 3D projection
   * \param problem :  the global problem to solve
   * \param localproblem :  the local problem to initialize
   */
  void fc3d_projection_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem);

  /** Initialize friction-contact 3D projection with regularization
   * \param problem :  the global problem to solve
   * \param localproblem :  the local problem to initialize
   */
  void fc3d_projection_initialize_with_regularization(FrictionContactProblem * problem, FrictionContactProblem * localproblem);

  /** Update friction-contact 3D projection solver: formalize local problem for one contact.
   * \param number (position in global matrix) of the considered contact
   * \param problem :  the global problem to solve
   * \param localproblem :  the local problem to initialize
   * \param reaction (only the block corresponding to the current contact will be modified,
   *   the rest is used to formalize the local problem)
   * \param options
   */
  void fc3d_projection_update(int number,  FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options);

  /** Update friction-contact 3D projection solver: formalize local problem for one contact.
   * \param number (position in global matrix) of the considered contact
   * \param problem :  the global problem to solve
   * \param localproblem :  the local problem to initialize
   * \param reaction (only the block corresponding to the current contact will be modified,
   * the rest is used to formalize the local problem)
   * \param options
   */
  void fc3d_projection_update_with_regularization(int number,  FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options);

  /** solve friction-contact 3D problem with projection assuming that M is diagonal
   * \param localproblem :  the local problem to initialize
   * \param reaction
   * \param options
   * \return 0 if successfull
   */
  int fc3d_projectionWithDiagonalization_solve(FrictionContactProblem * localproblem, double* reaction, SolverOptions * options);

  /** Update friction-contact 3D projection solver: formalize local problem for one contact.
   * \param number (position in global matrix) of the considered contact
   * \param problem :  the global problem to solve
   * \param localproblem :  the local problem to initialize
   * \param reaction (only the block corresponding to the current contact will be modified,
   * the rest is used to formalize the local problem)
   * \param options
   */
  void fc3d_projectionWithDiagonalization_update(int number,  FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options);

  int fc3d_projectionOnConeWithDiagonalization_setDefaultSolverOptions(SolverOptions* options);

  /** solve friction-contact 3D problem with projection on the Cone
   * \param localproblem :  the local problem to initialize
   * \param reaction
   * \param options
   * \return 0 if successfull
   */
  int fc3d_projectionOnCone_solve(FrictionContactProblem * localproblem, double* reaction, SolverOptions *options);

  int fc3d_projectionOnCone_setDefaultSolverOptions(SolverOptions* options);

  /** Update friction-contact 3D projection solver: formalize local problem for one contact.
   * \param number (position in global matrix) of the considered contact
   * \param problem :  the global problem to solve
   * \param localproblem :  the local problem to initialize
   * \param reaction (only the block corresponding to the current contact will be modified,
   * the rest is used to formalize the local problem)
   * \param options
  */
  void fc3d_projectionOnCylinder_update(int number,  FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options);

  int fc3d_projectionOnCylinder_setDefaultSolverOptions(SolverOptions* options);

  /** solve friction-contact 3D problem with projection on the Cone with local
   *   iteration up to convergence of the local problem
   * \param localproblem :  the local problem to initialize
   * \param reaction
   * \param options
   * \return 0 if successfull
   */
  int fc3d_projectionOnConeWithLocalIteration_solve(FrictionContactProblem * localproblem , double* reaction, SolverOptions * options);
  void fc3d_projectionOnConeWithLocalIteration_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options);
  void fc3d_projectionOnConeWithLocalIteration_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options);
  int fc3d_projectionOnConeWithLocalIteration_setDefaultSolverOptions(SolverOptions* options);
  /** solve friction-contact 3D problem with projection on the Cone
   * \param localproblem :  the local problem to initialize
   * \param reaction
   * \param options
   * \return 0 if successfull
   */
  int fc3d_projectionOnCone_velocity_solve(FrictionContactProblem * localproblem, double* reaction, SolverOptions * options);

  /** solve friction-contact 3D problem with projection on the (Tresca Cylinder)
   * \param localproblem :  the local problem to initialize
   * \param reaction
   * \param options
   * \return 0 if successfull
   */

  int fc3d_projectionOnCone_velocity_setDefaultSolverOptions(SolverOptions* options);
  
  int fc3d_projectionOnCylinder_solve(FrictionContactProblem * localproblem, double* reaction, SolverOptions * options);
  void fc3d_projectionOnCylinder_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* options);
  void fc3d_projectionOnCylinder_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options );

  /** solve friction-contact 3D problem with projection on the (Tresca Cylinder)
   * \param localproblem :  the local problem to initialize
   * \param reaction
   * \param options
   * \return 0 if successfull
   */
  int fc3d_projectionOnCylinderWithLocalIteration_solve(FrictionContactProblem * localproblem, double* reaction, SolverOptions * options);

  void fc3d_projectionOnCylinderWithLocalIteration_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* options, SolverOptions* localsolver_options);

  void fc3d_projectionOnCylinderWithLocalIteration_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options );
  
  int fc3d_projectionOnCylinderWithLocalIteration_setDefaultSolverOptions(SolverOptions* options);
  /** free memory for friction contact 3D projection solver
   * \param problem :  the  problem to free
   * \param localproblem :  the  problem to free
   * \param localsolver_options
   */
  void fc3d_projection_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options);

  /** free memory for friction contact 3D projection solver
   * \param problem :  the  problem to free
   * \param localproblem :  the  problem to free
   * \param localsolver_options
   */
  void fc3d_projection_with_regularization_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options);

  int fc3d_projectionOnConeWithRegularization_setDefaultSolverOptions(SolverOptions* options);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
