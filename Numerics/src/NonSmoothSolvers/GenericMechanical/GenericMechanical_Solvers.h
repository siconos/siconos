/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#ifndef GENERICMECHANICALSOLVERS_H
#define GENERICMECHANICALSOLVERS_H

/*!\file FrictionContact3D_Solvers.h
  \brief Subroutines for the resolution of contact problems with friction (3-dimensional case).\n

*/

/*! \page FC3DSolvers Friction-Contact 3D problems Solvers

This page gives an overview of the available solvers for friction-contact (3D) problems and their required parameters.

For each solver, the input argument are:
- a FrictionContactProblem
- the unknowns (reaction,velocity)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions structure, which handles iparam and dparam

\section fc3Dnsgs Non-Smooth Gauss Seidel Solver

\bf function: frictionContact3D_nsgs()
\bf parameters:


*/


#include "GenericMechanicalProblem.h"
#include "SolverOptions.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /** General interface to solvers for friction-contact 3D problem
  \param[in] , the structure which handles the generic mechanical problem
  \param[in-out] , reaction global vector (n)
  \param[in-out] , velocity global vector (n)
  \param[in,out] options structure used to define the solver(s) and their parameters
  \return result (0 if successful otherwise 1).
  */
  int genericMechanical_driver(GenericMechanicalProblem* problem, double *reaction , double *velocity, SolverOptions* options);

  GenericMechanicalProblem * buildEmptyGenericMechanicalProblem();
  void freeGenericMechanicalProblem(GenericMechanicalProblem * pGMP);
  /*return the localProblem (either lcp, linearSystem of fc3d*/
  void * addProblem(GenericMechanicalProblem * pGMP, int problemType, int size);
  void displayGMP(GenericMechanicalProblem * pGMP);
  void genericMechnicalProblem_setDefaultSolverOptions(SolverOptions* options, int id);
  void genericMechnical_printInFile(GenericMechanicalProblem*  problem, FILE* file);
  GenericMechanicalProblem* genericMechnical_newFromFile(FILE* file);
#ifdef __cplusplus
}
#endif

#endif
