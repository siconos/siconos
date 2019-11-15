/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
  \brief Subroutines for the resolution of contact problems.

*/

#include "GenericMechanicalProblem.h"
#include "SolverOptions.h"

#define NUMERICS_GMP_FREE_MATRIX 4
#define NUMERICS_GMP_FREE_GMP 8

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** General interface to solvers for friction-contact 3D problem
   * \param[in] problem the structure which handles the generic mechanical problem
   * \param[in,out]  reaction global vector (n)
   * \param[in,out]  velocity global vector (n)
   * \param[in,out] options structure used to define the solver(s) and their parameters
   *               option->iparam[0]:nb max of iterations
   *   option->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_LINESEARCH]:0 without 'LS' 1 with.
   *   option->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED]:0 GS block after block, 1 eliminate the equalities, 2 only one equality block, 3 solve the GMP as a MLCP.
   *   option->iparam[SICONOS_IPARAM_ITER_DONE]: output, number of GS it.
   *   options->dparam[SICONOS_DPARAM_TOL]: tolerance
   * \return result (0 if successful otherwise 1).
   */
  int gmp_driver(GenericMechanicalProblem* problem, double *reaction , double *velocity, SolverOptions* options);

  /** \addtogroup SetSolverOptions @{
   */
  void gmp_set_options(SolverOptions* options);
  /** @} */


  /* Alloc memory iff options->iWork and options->dWork  are  null.
   *\return  0 if the memory is not allocated. else return 1.
   */
  int gmp_working_memory_alloc(GenericMechanicalProblem* problem, SolverOptions* options);

  /* Free the Work memory, and set pointer to zero.*/
  void gmp_working_memory_free(GenericMechanicalProblem* problem, SolverOptions* options);

  /* Compute the error, return 0 iff success.*/
  int gmp_compute_error(GenericMechanicalProblem* pGMP, double *reaction , double *velocity, double tol, SolverOptions* options, double * err);

  /* Get the size of the double working zone memory.
   * \return the size of dWork
   */
  int gmp_get_nb_dwork(GenericMechanicalProblem* problem, SolverOptions* options);

  /* The  global Gauss-Seidel algorithm.
   */
  void gmp_gauss_seidel(GenericMechanicalProblem* pGMP, double * reaction, double * velocity, int * info, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
