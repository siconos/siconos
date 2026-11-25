/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef FRICTIONCONTACT3DPath_H
#define FRICTIONCONTACT3DPath_H

/*!\file fc3d_Path.h
  \brief Typedef and functions declarations related to NCP-Path solver for 3 dimension
  frictional contact problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for
  a specific contact
  - solve: solve the local problem
  - free

*/

#include "FrictionContactProblem.h"  // IWYU pragma: keep

#if defined(__cplusplus)
extern "C" {
#endif
/** Initialize friction-contact 3D Path solver
 * \param problem to solve
 * \param localproblem to solve
 * \param localsolver_options of the solver
 */
void fc3d_Path_initialize(FrictionContactProblem* problem,
                          FrictionContactProblem* localproblem,
                          SolverOptions* localsolver_options);

/** solve friction-contact 3D problem with Path
 * \param localproblem to solve
 * \param reaction
 * \param options of the solver
 * \return if sucessfull
 */
int fc3d_Path_solve(FrictionContactProblem* localproblem, double* reaction,
                    SolverOptions* options);

/** free memory for friction contact 3D Path solver */
void fc3d_Path_free(FrictionContactProblem*, SolverOptions*);

/**  compute error for  friction-contact 3D problem with Path
 * \param dimension of the global problem
 * \param[in,out] velocity vector (\warning in-out parameter )
 * \param reaction vector
 * \param output_error
 */
void fc3d_Path_computeError(int dimension, double* velocity, double* reaction,
                            double* output_error);

#if defined(__cplusplus)
}
#endif

#endif
