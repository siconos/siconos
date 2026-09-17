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

#ifndef GENERICMECHANICALSOLVERS_H
#define GENERICMECHANICALSOLVERS_H

/*!\file GenericMechanical_Solvers.h
  \brief Solvers for GenericMechanical problems.

  A GenericMechanical problem is a mixed non-smooth mechanical problem composed
  of several coupled sub-problems (equalities, Linear Complementarity Problems,
  Relay problems, 2D or 3D Friction-Contact problems, ...). The global matrix M
  is a block matrix whose diagonal blocks correspond to the local sub-problems.

  The main solver implemented here is a Non-Smooth Gauss-Seidel (NSGS) method:
  each block is solved sequentially while the values of the other blocks are
  fixed to their current iterate. Several variants are also provided for
  handling equality constraints (reduction, assembly, or MLCP reformulation).

  See \ref GenericMechanical_cst.h for the list of available solver ids and
  option parameters.

*/

#include "NumericsFwd.h"

#if defined(__cplusplus)
extern "C" {
#endif

/** General driver for GenericMechanical problems.
 *
 *  \param[in]     problem  the GenericMechanical problem to solve
 *  \param[in,out] reaction global reaction vector of size problem->globalSize
 *                          (input: initial guess, output: solution)
 *  \param[in,out] velocity global velocity vector of size problem->globalSize
 *                          (output: M*reaction + q)
 *  \param[in,out] options  solver options; must have been created for
 *                          SICONOS_GENERIC_MECHANICAL_NSGS and initialized with
 *                          gmp_set_default() or solver_options_create()
 *  \return 0 if successful, non-zero otherwise
 *
 *  The behavior is controlled by
 *  options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED]:
 *  - SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS: standard NSGS on all blocks
 *  - SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIES: eliminate equality blocks
 *    before NSGS
 *  - SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES: assemble all equality
 *    blocks into a single block before NSGS
 *  - SICONOS_GENERIC_MECHANICAL_MLCP_LIKE: reformulate as an MLCP (no FC3D)
 */
int gmp_driver(GenericMechanicalProblem* problem, double* reaction, double* velocity,
               SolverOptions* options);

/** \addtogroup SetSolverOptions
 * @{ */

/** Set default solver options for a GenericMechanical problem.
 *
 *  This fills \a options with default tolerances, iteration counts and default
 *  internal solvers for the four possible local problem types (LCP, FC3D,
 *  Relay, FC2D).
 *
 *  \param[in,out] options the SolverOptions structure to initialize
 */
void gmp_set_default(SolverOptions* options);

/** @} */

/** Allocate working memory if not already present.
 *
 *  This function allocates options->dWork only if it is currently NULL.
 *
 *  \param[in]     problem the GenericMechanical problem
 *  \param[in,out] options the solver options
 *  \return 0 if memory was already allocated, 1 if it has been allocated
 */
int gmp_working_memory_alloc(GenericMechanicalProblem* problem, SolverOptions* options);

/** Compute the global residual of a GenericMechanical problem.
 *
 *  The local problems are updated with the current \a reaction iterate and the
 *  corresponding local errors are evaluated. The maximum local error (scaled
 *  appropriately) is returned in \a err.
 *
 *  \param[in]  problem  the GenericMechanical problem
 *  \param[in]  reaction current global reaction iterate
 *  \param[out] velocity global velocity vector, M*reaction + q
 *  \param[in]  tol      requested tolerance
 *  \param[in]  options  solver options
 *  \param[out] err      computed maximum error
 *  \return 0 if err <= tol, 1 otherwise
 */
int gmp_compute_error(const GenericMechanicalProblem* problem, double* reaction,
                      double* velocity, double tol, SolverOptions* options, double* err);

/** Return the size (number of doubles) needed for options->dWork.
 *
 *  \param[in] problem the GenericMechanical problem
 *  \param[in] options the solver options
 *  \return the required size of the double working array
 */
int gmp_get_nb_dwork(GenericMechanicalProblem* problem, SolverOptions* options);

/** Non-Smooth Gauss-Seidel iteration for a GenericMechanical problem.
 *
 *  Sequentially solve each local sub-problem while freezing the other blocks.
 *  The local right-hand side is built from the global q and the off-diagonal
 *  block row product M_offdiag * reaction.
 *
 *  \param[in]     problem      the GenericMechanical problem
 *  \param[in,out] reaction  global reaction vector (initial guess / solution)
 *  \param[in,out] velocity  global velocity vector
 *  \param[out]    info      0 if convergence, non-zero otherwise
 *  \param[in,out] options   solver options
 */
void gmp_gauss_seidel(GenericMechanicalProblem* problem, double* reaction, double* velocity,
                      int* info, SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
