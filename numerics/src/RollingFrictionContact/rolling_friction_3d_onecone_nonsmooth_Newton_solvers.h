/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef ROLLING_FRICTION_onecone_nonsmooth_Newton_solvers_H
#define ROLLING_FRICTION_onecone_nonsmooth_Newton_solvers_H

/*!\file rolling_friction_3d_onecone_nonsmooth_Newton_solvers.h
  \brief Typedef and functions declarations related to Newton solver for rolling friction 3D
  problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for
  a specific contact
  - solve: solve the local problem
  - free

*/

#include "NumericsFwd.h"    // for Plasticity2DProblem, SolverOptions
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

typedef void (*computeNonsmoothFunction)(double*, double*, double, double, double*, double*,
                                         double*, double*);

/** initialize rolling friction 3D Newton solver
 * \param problem to solve
 * \param localproblem to solve
 * \param options of the solver
 */
void rolling_friction_3d_onecone_nonsmooth_Newton_solvers_initialize(RollingFrictionContactProblem* problem,
                                                      RollingFrictionContactProblem* localproblem,
                                                      SolverOptions* options);

/** solve rolling friction 3D problem with Newton
 * \param localproblem to solve
 * \param options of the solver
 * \return 0 iff successful.
 */
int rolling_friction_3d_onecone_nonsmooth_Newton_solvers_solve(RollingFrictionContactProblem* localproblem, double*,
                                                SolverOptions* options);

/** free memory for rolling friction 3D Newton solver
    \param problem the global problem to solve
    \param localproblem for freeing matrix0
    \param localsolver_options options of the solver
 */
void rolling_friction_3d_onecone_nonsmooth_Newton_solvers_free(RollingFrictionContactProblem* problem,
                                                RollingFrictionContactProblem* localproblem,
                                                SolverOptions* localsolver_options);

/** compute error for rolling friction 3D problem with Newton
 *  \param dimension of the global problem
 *  \param[in,out] velocity vector
 *  \param reaction global reaction vector
 *  \param output_error
 */
void rolling_friction_3d_onecone_nonsmooth_Newton_solvers_computeError(int dimension, double* velocity,
                                                        double* reaction,
                                                        double* output_error);

/** Update rolling friction 3D problem: formalize local problem for one contact
    \param problem the global problem to solve
    \param localproblem the local problem to solve
    \param number (position in global matrix) of the considered contact
    \param reaction global reaction (only the block corresponding to the
    current contact will be modified
    \param options of the solver

    the rest is used to formalize the local problem)
*/
void rolling_friction_3d_onecone_nonsmooth_Newton_update(int number, RollingFrictionContactProblem* problem,
                                             RollingFrictionContactProblem* localproblem,
                                             double* reaction, SolverOptions* options);

int rolling_friction_3d_onecone_nonsmooth_Newton_solvers_solve_direct(RollingFrictionContactProblem* localproblem,
                                                       double* R, SolverOptions* options);

int rolling_friction_3d_onecone_nonsmooth_Newton_solvers_solve_damped(RollingFrictionContactProblem* localproblem,
                                                       double* R, SolverOptions* options);

int rolling_friction_3d_onecone_nonsmooth_Newton_solvers_solve_hybrid(RollingFrictionContactProblem* localproblem,
                                                       double* local_reaction,
                                                       SolverOptions* options);

/** Set default options for onecone NSN solver
 * \param options the solver options
 */
void rolling_friction_3d_onecone_nsn_set_default(SolverOptions* options);

/** Set default options for onecone NSN GP solver
 * \param options the solver options
 */
void rolling_friction_3d_onecone_nsn_gp_set_default(SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
