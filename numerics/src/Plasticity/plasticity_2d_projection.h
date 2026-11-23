// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef PLASTICITY2DProjection_H
#define PLASTICITY2DProjection_H

/*!\file plasticity_2d_projection.h
  \brief Typedef and functions declarations related to projection solver for
  Mohr Coulomb 2D

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
#include "NumericsFwd.h"  // for PlasticityProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

/** Initialize Mohr Coulomb 2D projection
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 */
void plasticity_2d_projection_initialize(PlasticityProblem* problem,
                                         PlasticityProblem* localproblem);

/** Initialize Mohr Coulomb 2D projection with regularization
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 */
void plasticity_2d_projection_initialize_with_regularization(PlasticityProblem* problem,
                                                             PlasticityProblem* localproblem);

/** Update Mohr Coulomb 2D projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param reaction (only the block corresponding to the current contact will be modified,
 *   the rest is used to formalize the local problem)
 * \param options
 */
void plasticity_2d_projection_update(int number, PlasticityProblem* problem,
                                     PlasticityProblem* localproblem, double* reaction,
                                     SolverOptions* options);

/** solve Mohr Coulomb 2D problem with projection on the Cone
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int plasticity_2d_projectionOnCone_solve(PlasticityProblem* localproblem, double* reaction,
                                         SolverOptions* options);

/** solve Mohr Coulomb 2D problem with projection on the Cone with local
 *   iteration up to convergence of the local problem
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int plasticity_2d_projectionOnConeWithLocalIteration_solve(PlasticityProblem* localproblem,
                                                           double* reaction,
                                                           SolverOptions* options);
void plasticity_2d_projectionOnConeWithLocalIteration_free(PlasticityProblem* problem,
                                                           PlasticityProblem* localproblem,
                                                           SolverOptions* localsolver_options);
void plasticity_2d_projectionOnConeWithLocalIteration_initialize(
    PlasticityProblem* problem, PlasticityProblem* localproblem,
    SolverOptions* localsolver_options);
/** free memory for Mohr Coulomb 2D projection solver
 * \param problem :  the  problem to free
 * \param localproblem :  the  problem to free
 * \param localsolver_options
 */
void plasticity_2d_projection_free(PlasticityProblem* problem, PlasticityProblem* localproblem,
                                   SolverOptions* localsolver_options);

/* ============================================================================
 * Von Mises Projection Solver
 * ============================================================================ */

/** Von Mises (J2) one-element problem (cone coordinates, see
 * plasticity_2d_vonmises_projection.h): one projection step
 *   sigma <- P_C(sigma - rho (M sigma + q)),  rho = 1 / M_NN.
 * \param localproblem the local problem (one element, model Von Mises)
 * \param stress_local the local stress (input/output)
 * \param options solver options (unused)
 * \return 0
 */
int plasticity_2d_projectionOnVonMises_solve(PlasticityProblem* localproblem,
                                             double* stress_local, SolverOptions* options);

/** Von Mises (J2) one-element problem: exact solution
 * (plasticity_2d_vonMises_local_solve).
 * \param localproblem the local problem (one element, model Von Mises)
 * \param stress_local the local stress (output)
 * \param options solver options (unused)
 * \return 0 if successful, 1 if the local matrix is not suitable (M_NN <= 0, S singular)
 */
int plasticity_2d_vonMises_exact_solve(PlasticityProblem* localproblem, double* stress_local,
                                       SolverOptions* options);

/** Set default solver options for Von Mises projection
 * \param options the solver options to set
 */
void plasticity_2d_vonmises_set_default(SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
