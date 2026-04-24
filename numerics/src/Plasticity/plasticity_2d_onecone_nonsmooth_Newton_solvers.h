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
#ifndef PLASTICITY2D_onecone_nonsmooth_Newton_solvers_H
#define PLASTICITY2D_onecone_nonsmooth_Newton_solvers_H

/*!\file plasticity_2d_onecone_nonsmooth_Newton_solvers.h
  \brief Typedef and functions declarations related to Newton solver for Mohr Coulomb 2D
  problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for
  a specific contact
  - solve: solve the local problem
  - free

*/

#include "NumericsFwd.h"  // for PlasticityProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

typedef void (*computeNonsmoothFunction)(double*, double*, double, double, double*, double*,
                                         double*, double*);

/** initialize Mohr Coulomb 2D Newton solver
 * \param problem to solve
 * \param localproblem to solve
 * \param options of the solver
 */
void plasticity_2d_onecone_nonsmooth_Newton_solvers_initialize(PlasticityProblem* problem,
                                                               PlasticityProblem* localproblem,
                                                               SolverOptions* options);

/** solve Mohr Coulomb 2D problem with Newton
 * \param localproblem to solve
 * \param options of the solver
 * \return 0 iff successful.
 */
int plasticity_2d_onecone_nonsmooth_Newton_solvers_solve(PlasticityProblem* localproblem,
                                                         double*, SolverOptions* options);

/** free memory for Mohr Coulomb 2D Newton solver
    \param problem the global problem to solve
    \param localproblem for freeing matrix0
    \param localsolver_options options of the solver
 */
void plasticity_2d_onecone_nonsmooth_Newton_solvers_free(PlasticityProblem* problem,
                                                         PlasticityProblem* localproblem,
                                                         SolverOptions* localsolver_options);

/** compute error for Mohr Coulomb 2D problem with Newton
 *  \param dimension of the global problem
 *  \param[in,out] velocity vector
 *  \param reaction global reaction vector
 *  \param output_error
 */
void plasticity_2d_onecone_nonsmooth_Newton_solvers_computeError(int dimension,
                                                                 double* velocity,
                                                                 double* reaction,
                                                                 double* output_error);

/** Update Mohr Coulomb 2D problem: formalize local problem for one contact
    \param problem the global problem to solve
    \param localproblem the local problem to solve
    \param number (position in global matrix) of the considered contact
    \param reaction global reaction (only the block corresponding to the
    current contact will be modified
    \param options of the solver

    the rest is used to formalize the local problem)
*/
void plasticity_2d_onecone_nonsmooth_Newton_AC_update(int number, PlasticityProblem* problem,
                                                      PlasticityProblem* localproblem,
                                                      double* reaction,
                                                      SolverOptions* options);

int plasticity_2d_onecone_nonsmooth_Newton_solvers_solve_direct(
    PlasticityProblem* localproblem, double* R, SolverOptions* options);

int plasticity_2d_onecone_nonsmooth_Newton_solvers_solve_damped(
    PlasticityProblem* localproblem, double* R, SolverOptions* options);

int plasticity_2d_onecone_nonsmooth_Newton_solvers_solve_hybrid(
    PlasticityProblem* localproblem, double* local_reaction, SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
