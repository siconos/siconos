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
#ifndef MOHRCOULOMB2DProjection_H
#define MOHRCOULOMB2DProjection_H

/*!\file mc2d_projection.h
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
#include "NumericsFwd.h"    // for MohrCoulomb2DProblem, SolverOptions
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/** Initialize Mohr Coulomb 2D projection
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 */
void mc2d_projection_initialize(MohrCoulomb2DProblem* problem,
                                MohrCoulomb2DProblem* localproblem);

/** Initialize Mohr Coulomb 2D projection with regularization
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 */
void mc2d_projection_initialize_with_regularization(MohrCoulomb2DProblem* problem,
                                                    MohrCoulomb2DProblem* localproblem);

/** Update Mohr Coulomb 2D projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param reaction (only the block corresponding to the current contact will be modified,
 *   the rest is used to formalize the local problem)
 * \param options
 */
void mc2d_projection_update(int number, MohrCoulomb2DProblem* problem,
                            MohrCoulomb2DProblem* localproblem, double* reaction,
                            SolverOptions* options);

/** Update Mohr Coulomb 2D projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param reaction (only the block corresponding to the current contact will be modified,
 * the rest is used to formalize the local problem)
 * \param options
 */
void mc2d_projection_update_with_regularization(int number, MohrCoulomb2DProblem* problem,
                                                MohrCoulomb2DProblem* localproblem,
                                                double* reaction, SolverOptions* options);

/** solve Mohr Coulomb 2D problem with projection assuming that M is diagonal
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int mc2d_projectionWithDiagonalization_solve(MohrCoulomb2DProblem* localproblem,
                                             double* reaction, SolverOptions* options);

/** Update Mohr Coulomb 2D projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param reaction (only the block corresponding to the current contact will be modified,
 * the rest is used to formalize the local problem)
 * \param options
 */
void mc2d_projectionWithDiagonalization_update(int number, MohrCoulomb2DProblem* problem,
                                               MohrCoulomb2DProblem* localproblem,
                                               double* reaction, SolverOptions* options);

/** solve Mohr Coulomb 2D problem with projection on the Cone
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int mc2d_projectionOnCone_solve(MohrCoulomb2DProblem* localproblem, double* reaction,
                                SolverOptions* options);

/** Update Mohr Coulomb 2D projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param reaction (only the block corresponding to the current contact will be modified,
 * the rest is used to formalize the local problem)
 * \param options
 */
void mc2d_projectionOnCylinder_update(int number, MohrCoulomb2DProblem* problem,
                                      MohrCoulomb2DProblem* localproblem, double* reaction,
                                      SolverOptions* options);

/** solve Mohr Coulomb 2D problem with projection on the Cone with local
 *   iteration up to convergence of the local problem
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int mc2d_projectionOnConeWithLocalIteration_solve(MohrCoulomb2DProblem* localproblem,
                                                  double* reaction, SolverOptions* options);
void mc2d_projectionOnConeWithLocalIteration_free(MohrCoulomb2DProblem* problem,
                                                  MohrCoulomb2DProblem* localproblem,
                                                  SolverOptions* localsolver_options);
void mc2d_projectionOnConeWithLocalIteration_initialize(MohrCoulomb2DProblem* problem,
                                                        MohrCoulomb2DProblem* localproblem,
                                                        SolverOptions* localsolver_options);
/** free memory for Mohr Coulomb 2D projection solver
 * \param problem :  the  problem to free
 * \param localproblem :  the  problem to free
 * \param localsolver_options
 */
void mc2d_projection_free(MohrCoulomb2DProblem* problem,
                          MohrCoulomb2DProblem* localproblem,
                          SolverOptions* localsolver_options);

/** free memory for Mohr Coulomb 2D projection solver
 * \param problem :  the  problem to free
 * \param localproblem :  the  problem to free
 * \param localsolver_options
 */
void mc2d_projection_with_regularization_free(MohrCoulomb2DProblem* problem,
                                              MohrCoulomb2DProblem* localproblem,
                                              SolverOptions* localsolver_options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
