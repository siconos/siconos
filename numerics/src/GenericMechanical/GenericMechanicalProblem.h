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

/*! \file GenericMechanicalProblem.h
 * \brief struct for GenericMechanicalProblem
 */

#ifndef NUMERICSGENERICMECHANICALPROBLEM_H
#define NUMERICSGENERICMECHANICALPROBLEM_H

#include "NumericsFwd.h"
#include <stdio.h>
/* void * solverFC3D; */
/* void * solverEquality; */
/* void * solverLCP; */
/* void * solverMLCP; */

/** \struct GenericMechanicalProblem GenericMechanicalProblem.h
 *  \param numberOfBlockLine The number of  line of blocks.
 *   \param M a sparse blocks matrix.
 *   \param q a dense vector.
 *   \param size sizes of the local problems (needed in the dense case)
 *   \param nextProblem the list of the next problems
 *   \param prevProblem the list of the previous problems
 *   Remark:
 *   The M and q contains the matrices of the GMP problem. The sub problems (problems) has also a M and q member usfull for the computation of the local error.
 *
 */
struct listNumericsProblem
{
  int type;
  void * problem;
  double *q;/*a pointer on the q of the problem*/
  int size;/*size of the local problem.(needed because of dense case)*/
  int error;/*non-zero if there was an error reported*/
  struct listNumericsProblem * nextProblem;
  struct listNumericsProblem * prevProblem;
};


/** \struct GenericMechanicalProblem GenericMechanicalProblem.h
 * \param numberOfBlockLine The number of  line of blocks.
 * \param M : NumericsMatrix sparseblock matrix set by the user
 * \param q : dense vector set by the user
 * \param size : maximal size of local problem
 * \param maxLocalSize "private" manage by gmp_add
 * \param firstListElem "private" manage by gmp_add
 * \param lastListElem  "private" manage by gmp_add
 *
 *  Remark:
 *  The M and q contains the matrices of the GMP problem.
 *  The sub problems (problems) has also a M and q member usfull for the computation of the local error.
 *
 * ONLY q and M must be allocated/free by the users, the others fields are private:
 * DO NOT FILL THIS STRUCTURE BY YOURSELF, BUT USE THE
 * - genericMechanicalProblem_new() ,
 * - gmp_add() ,
 * - and genericMechanicalProblem_free() FUNCTIONS.
 */
struct GenericMechanicalProblem
{
  /*Number of line of blocks.*/
  /*PRIVATE: manage by gmp_add.*/
  int size;
  /*maximal size of local problem.*/
  /*PRIVATE: manage by gmp_add.*/
  int maxLocalSize;
  /*must be set by the user.*/
  NumericsMatrix* M;
  /*must be set by the user.*/
  double* q;
  /*PRIVATE: manage by gmp_add.*/
  listNumericsProblem *firstListElem;
  /*PRIVATE: manage by gmp_add.*/
  listNumericsProblem *lastListElem;
  //  void * * problems;
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /* Build an empty GenericMechanicalProblem
   * \return a pointer on the built GenericMechanicalProblem.
   */
  GenericMechanicalProblem * genericMechanicalProblem_new(void);

  /* Free the list of the contained sub-problem, coherently with the memory allocated in the gmp_add function, it also free the pGMP.
   */
  void genericMechanicalProblem_free(GenericMechanicalProblem * pGMP, unsigned int level);

  /* To print a GenericMechanicalProblem in a file.
   *  \param[in] problem, the printed problem.
   *  \param[in,out] output file.
   */
  void genericMechanicalProblem_printInFile(GenericMechanicalProblem*  problem, FILE* file);

  /*To build a GenericMechanicalProblem from a file.
   * \parm[in] file, a file containing the GenericMechanicalProblem.
   * \return the GenericMechanicalProblem.
   */
  GenericMechanicalProblem* genericMechanical_newFromFile(FILE* file);

  /* A recursive displaying method.
   *  \param[in], pGMP the displayed problem.
   */
  void genericMechanicalProblem_display(GenericMechanicalProblem * pGMP);

  /* Insert a problem in the GenericMechanicalProblem pGMP. The memory of the elematary block is not managed. The user has to ensure it.
   *   In the case of SICONOS, the Kernel ensure this allocation in building the global problem. In other words, the matrix0 is shared with the global NumericsMatrix,
   * the plug is done in the function gmp_gauss_seidel (ie: localProblem->M->matrix0= m->block[diagBlockNumber];)
   * \param[in,out] pGMP a pointer.
   * \param[in] problemType type of the added sub-problem (either SICONOS_NUMERICS_PROBLEM_LCP, SICONOS_NUMERICS_PROBLEM_EQUALITY, SICONOS_NUMERICS_PROBLEM_FC3D, or SICONOS_NUMERICS_PROBLEM_RELAY)
   * \param[in] size size of the formulation (dim of the LCP, or dim of the linear system, 3 for the fc3d)
   * \ return the localProblem (either lcp, linearSystem of fc3d
   */
  void * gmp_add(GenericMechanicalProblem * pGMP, int problemType, int size);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
