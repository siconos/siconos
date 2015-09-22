/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef SOCLCPProjection_H
#define SOCLCPProjection_H

/*!\file soclcp_projection.h
  \brief Typedef and functions declarations related to projection solver for SOCLCP

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
#include "SOCLCP_Solvers.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

/** Initialize SOCLCP projection
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param options 
 */
void soclcp_projection_initialize(SecondOrderConeLinearComplementarityProblem * problem,
                                  SecondOrderConeLinearComplementarityProblem * localproblem,
                                  SolverOptions * options);

/** Initialize SOCLCP projection with regularization
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 */
void soclcp_projection_initialize_with_regularization(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem);

/** Update SOCLCP Projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param r (only the block corresponding to the current contact will be modified,
 *   the rest is used to formalize the local problem)
 * \param options
 */
void soclcp_projection_update(int number,  SecondOrderConeLinearComplementarityProblem* problem, SecondOrderConeLinearComplementarityProblem* localproblem, double* r, SolverOptions* options);

/** Update SOCLCP projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param r (only the block corresponding to the current contact will be modified,
 * the rest is used to formalize the local problem)
 * \param options
 */
void soclcp_projection_update_with_regularization(int number,  SecondOrderConeLinearComplementarityProblem* problem, SecondOrderConeLinearComplementarityProblem* localproblem, double* r, SolverOptions* options);

/** solve SOCLCP problem with projection on the Cone
 * \param localproblem :  the local problem to initialize
 * \param r
 * \param options
 * \return 0 if successfull
 */
int soclcp_projectionOnCone_solve(SecondOrderConeLinearComplementarityProblem * localproblem, double* r, SolverOptions *options);

/** solve SOCLCP problem with projection on the Cone with local
 *   iteration up to convergence of the local problem
 * \param localproblem :  the local problem to initialize
 * \param r
 * \param options
 * \return 0 if successfull
 */
int soclcp_projectionOnConeWithLocalIteration_solve(SecondOrderConeLinearComplementarityProblem * localproblem , double* r, SolverOptions * options);
void soclcp_projectionOnConeWithLocalIteration_free(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, SolverOptions* localsolver_options);
void soclcp_projectionOnConeWithLocalIteration_initialize(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, SolverOptions* localsolver_options);

/** solve SOCLCP problem with projection on the Cone
 * \param localproblem :  the local problem to initialize
 * \param r
 * \param options
 * \return 0 if successfull
 */
int soclcp_projectionOnCone_v_solve(SecondOrderConeLinearComplementarityProblem * localproblem, double* r, SolverOptions * options);

/** solve SOCLCP problem with projection on the (Tresca Cylinder)
 * \param localproblem :  the local problem to initialize
 * \param r
 * \param options
 * \return 0 if successfull
 */
int soclcp_projectionOnCylinder_solve(SecondOrderConeLinearComplementarityProblem * localproblem, double* r, SolverOptions * options);

/** free memory for friction contact 3D projection solver
 * \param problem :  the  problem to free
 * \param localproblem :  the  problem to free
 * \param localsolver_options
 */
void soclcp_projection_free(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, SolverOptions* localsolver_options);

/** free memory for friction contact 3D projection solver
 * \param problem :  the  problem to free
 * \param localproblem :  the  problem to free
 * \param localsolver_options
 */
void soclcp_projection_with_regularization_free(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, SolverOptions* localsolver_options);

int soclcp_projection_setDefaultSolverOptions(SolverOptions* options);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
