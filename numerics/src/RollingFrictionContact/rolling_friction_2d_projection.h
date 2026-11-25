// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef ROLLINGFRICTIONCONTACT2DProjection_H
#define ROLLINGFRICTIONCONTACT2DProjection_H

/*!\file rolling_friction_2d_projection.h
  \brief Typedef and functions declarations related to projection solver for 3 dimension
  frictional contact problems

  Each solver must have 4 functions in its interface:
  - initialize: link global static variables to the considered problem (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for
  a specific contact
  - solve: solve the local problem
  - free
  We consider a "global" (ie for several contacts) problem, used to initialize the static
  global variables. Then a "local" (ie for one contact => size = 3) problem is built (update
  function) and solved (solve function).

  Two different storages are available for M: dense and sparse block.

*/

#include "NumericsFwd.h"  // for RollingFrictionContactProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif
/** Update friction-contact 3D projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param reaction (only the block corresponding to the current contact will be modified,
 *   the rest is used to formalize the local problem)
 * \param options
 */
void rolling_friction_2d_projection_update(int number, RollingFrictionContactProblem* problem,
                                           RollingFrictionContactProblem* localproblem,
                                           double* reaction, SolverOptions* options);

void rolling_friction_2d_projection_initialize(RollingFrictionContactProblem* problem,
                                               RollingFrictionContactProblem* localproblem);
void rolling_friction_2d_projection_free(RollingFrictionContactProblem* problem,
                                         RollingFrictionContactProblem* localproblem,
                                         SolverOptions* localsolver_options);

int rolling_friction_2d_projectionOnCone_solve(RollingFrictionContactProblem* localproblem,
                                               double* reaction, SolverOptions* options);
int rolling_friction_2d_projectionOnCone_setDefaultSolverOptions(SolverOptions* options);

/** solve friction-contact 3D problem with projection on the Cone with local
 *   iteration up to convergence of the local problem
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int rolling_friction_2d_projectionOnConeWithLocalIteration_solve(
    RollingFrictionContactProblem* localproblem, double* reaction, SolverOptions* options);
void rolling_friction_2d_projectionOnConeWithLocalIteration_free(
    RollingFrictionContactProblem* problem, RollingFrictionContactProblem* localproblem,
    SolverOptions* localsolver_options);
void rolling_friction_2d_projectionOnConeWithLocalIteration_initialize(
    RollingFrictionContactProblem* problem, RollingFrictionContactProblem* localproblem,
    SolverOptions* localsolver_options);
int rolling_friction_2d_projectionOnConeWithLocalIteration_setDefaultSolverOptions(
    SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
