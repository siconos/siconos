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

extern unsigned int NUMERICS_GMP_FREE_MATRIX;
extern unsigned int NUMERICS_GMP_FREE_GMP;

#ifdef __cplusplus
extern "C"
{
#endif

  /** General interface to solvers for friction-contact 3D problem
  \param[in] , the structure which handles the generic mechanical problem
  \param[in-out] , reaction global vector (n)
  \param[in-out] , velocity global vector (n)
  \param[in,out] options structure used to define the solver(s) and their parameters
                 option->iparam[0]:nb max of iterations
     option->iparam[1]:0 without 'LS' 1 with.
     option->iparam[2]:0 GS block after block, 1 eliminate the equalities, 2 only one equality block.
     option->iparam[3]: output, number of GS it.
     options->dparam[0]: tolerance
  \return result (0 if successful otherwise 1).
  */
  int genericMechanical_driver(GenericMechanicalProblem* problem, double *reaction , double *velocity, SolverOptions* options);
  /* Build an empty GenericMechanicalProblem
     \return a pointer on the built GenericMechanicalProblem.
   */
  GenericMechanicalProblem * buildEmptyGenericMechanicalProblem();
  /* Free the list of the contained sub-problem, coherently with the memory allocated in the addProblem function, it also free the pGMP.
   */
  void freeGenericMechanicalProblem(GenericMechanicalProblem * pGMP, unsigned int level);
  /* Insert a problem in the GenericMechanicalProblem pGMP. The memory of the elematary block is not managed. The user has to ensure it.
     In the case of SICONOS, the Kernel ensure this allocation in building the global problem. In other words, the matrix0 is shared with the global NumericsMatrix,
     the plug is done in the function genericMechanicalProblem_GS (ie: localProblem->M->matrix0= m->block[diagBlockNumber];)
   \param[in-out], pGMP a pointer.
   \param[in], problemType, type of the added sub-problem (either SICONOS_NUMERICS_PROBLEM_LCP, SICONOS_NUMERICS_PROBLEM_EQUALITY or SICONOS_NUMERICS_PROBLEM_FC3D)
   \param[in], size of the formulation (dim of the LCP, or dim of the linear system, 3 for the fc3d)
   \ return the localProblem (either lcp, linearSystem of fc3d
   */
  void * addProblem(GenericMechanicalProblem * pGMP, int problemType, int size);
  /*A recursive displaying method.
    \param[in], pGMP the displayed problem.
   */
  void displayGMP(GenericMechanicalProblem * pGMP);
  /*Builder of options, it supposes that options is not NULL.
   \param [in] id, not used in current version
   \param [in-out] options, the filled options.(by default LEMKE and Quartic)
   */
  void genericMechnicalProblem_setDefaultSolverOptions(SolverOptions* options, int id);
  /*To print a GenericMechanicalProblem in a file.
    \param[in] problem, the printed problem.
    \param[in-out] output file.
   */
  void genericMechnical_printInFile(GenericMechanicalProblem*  problem, FILE* file);
  /*To build a GenericMechanicalProblem from a file.
    \parm[in] file, a file containing the GenericMechanicalProblem.
    \return the built GenericMechanicalProblem.
   */
  GenericMechanicalProblem* genericMechnical_newFromFile(FILE* file);

  /*Alloc memory iff options->iWork options->dWork and are  null.
  Return 0 if the memory is not allocated. else return 1.*/
  int genericMechnical_alloc_working_memory(GenericMechanicalProblem* problem, SolverOptions* options);
  /*free the Work memory, and set pointer to zero.*/
  void genericMechnical_free_working_memory(GenericMechanicalProblem* problem, SolverOptions* options);

  int genericMechnical_getNbDWork(GenericMechanicalProblem* problem, SolverOptions* options);
  void genericMechanicalProblem_GS(GenericMechanicalProblem* pGMP, double * reaction, double * velocity, int * info, SolverOptions* options);
#ifdef __cplusplus
}
#endif

#endif
