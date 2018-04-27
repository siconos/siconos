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
 * \param maxLocalSize "private" manage by addProblem
 * \param firstListElem "private" manage by addProblem
 * \param lastListElem  "private" manage by addProblem
 *
 *  Remark:
 *  The M and q contains the matrices of the GMP problem.
 *  The sub problems (problems) has also a M and q member usfull for the computation of the local error.
 *
 * ONLY q and M must be allocated/free by the users, the others fields are private:
 * DO NOT FILL THIS STRUCTURE BY YOURSELF, BUT USE THE
 * - buildEmptyGenericMechanicalProblem() ,
 * - addProblem() ,
 * - and freeGenericMechanicalProblem() FUNCTIONS.
 */
struct GenericMechanicalProblem
{
  /*Number of line of blocks.*/
  /*PRIVATE: manage by addProblem.*/
  int size;
  /*maximal size of local problem.*/
  /*PRIVATE: manage by addProblem.*/
  int maxLocalSize;
  /*must be set by the user.*/
  NumericsMatrix* M;
  /*must be set by the user.*/
  double* q;
  /*PRIVATE: manage by addProblem.*/
  listNumericsProblem *firstListElem;
  /*PRIVATE: manage by addProblem.*/
  listNumericsProblem *lastListElem;
  //  void * * problems;
};




#endif
