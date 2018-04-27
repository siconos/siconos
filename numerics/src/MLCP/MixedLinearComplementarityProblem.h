/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#ifndef MLCP_PROBLEM_H
#define MLCP_PROBLEM_H

/*!\file MixedLinearComplementarityProblem.h
  \brief Structure used to define a Mixed Linear Complementarity Problem

  \author Olivier Bonnefon
*/

#include <stdio.h>
#include "NumericsFwd.h"
#include "SiconosConfig.h"

/** \struct MixedLinearComplementarityProblem MixedLinearComplementarityProblem.h
 *  The Structure that contains and defines  \ref MLCProblem. Find \f$(z,w)\f$ such that:\n
 * \f$
 * \left\{
 * \begin{array}{l}
 *  M \ z + q = w \\
 * w_1=0 \\
 * 0 \le w_{2} \perp v \ge 0
 * \end{array}
 * \right.
 * \text{ with }
 * z=
 * \left[
 * \begin{array}{c}
 * u\\
 * v\\
 * \end{array}
 * \right]
 * \text{ and }
 * w=
 * \left[
 * \begin{array}{c}
 * w_{1}\\
 * w_{2}\\
 * \end{array}
 * \right]
 * \f$
 * \f$ u, w_{1}\f$ are vectors of size n.\n
 * \f$ v, w_{2}\f$ are vectors of size m.
 * See \ref MLCProblem for more details.
 */
struct MixedLinearComplementarityProblem
{
  int isStorageType1; /**< boolean for storageType1 1 if the problem
                         is saved using (M,q),  0 otherwise */
  int isStorageType2; /**< boolean for storageType2 1 if the problem
                         is saved using (A,B,C,D,a,b), 0 otherwise*/
  int n; /**< number of equality constraints */
  int m; /**< number of complementarity constraints */
  int * blocksRows;  /**< The rows from blocksRows[i] to blocksRows[i+1]-1
                        forms a block of equalities iif bloksIsComp[i]=0,
                        else the block is a complementarity block.
                        The number of total blocks is given by NbBlocks
                        such that blocksRows[NbBlocks] = n+m */
  int * blocksIsComp; /**< if bloksIsComp[i]=0, then block i formed by the rows
                         from blocksRows[i] to blocksRows[i+1]-1 is an equality block
                         else the block is a complementarity block.
                      */
  NumericsMatrix* M; /**< M matrix of the MLCP */
  double *q; /**< q vector of the MLCP */
  /** NumericsMatrix* Bblock;*/ /**< Bblock  ?*/
  double *A; /**< A matrix of the MLCP */
  double *B; /**< B matrix of the MLCP */
  double *C; /**< C matrix of the MLCP */
  double *D; /**< D matrix of the MLCP */
  double *a; /**< a vector of the MLCP */
  double *b; /**< b vector of the MLCP */
};


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void mixedLinearComplementarity_display(MixedLinearComplementarityProblem* p);

  /** \fn int mixedLinearComplementarity_printInFile(MixedLinearComplementarityProblem*  problem, FILE* file)
   *  \brief function to write in a file a MixedLinearComplementarityProblem
   *  \param problem pointer to a MixedLinearComplementarityProblem to print
   *  \param file pointer to a FILE
   *  \return 0 if ok
   */
  int mixedLinearComplementarity_printInFile(MixedLinearComplementarityProblem*  problem, FILE* file);


  /** \fn  int mixedLinearComplementarity_newFromFile(MixedLinearComplementarityProblem* problem, FILE* file)
   *  \brief function to read and create a MixedLinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a MixedLinearComplementarityProblem to create
   *  \param file pointer to a FILE
   *  \return 0 if ok
   */
  int mixedLinearComplementarity_newFromFile(MixedLinearComplementarityProblem* problem, FILE* file);

  /** \fn  int mixedLinearComplementarity_newFromFileOld(MixedLinearComplementarityProblem* problem, FILE* file)
   *  \brief function to read and create a MixedLinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a MixedLinearComplementarityProblem to create
   *  \param file pointer to a FILE
   *  \return 0 if ok
   */
  int mixedLinearComplementarity_newFromFileOld(MixedLinearComplementarityProblem* problem, FILE* file);

  /** \fn  int mixedLinearComplementarity_newFromFilename(MixedLinearComplementarityProblem* problem, FILE* MLCPfile)
   *  \brief function to read and create a MixedLinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a MixedLinearComplementarityProblem to create
   *  \param filename that contains the mlcp
   *  \return 0 if ok
   */
  int mixedLinearComplementarity_newFromFilename(MixedLinearComplementarityProblem* problem, char* filename);

  /** \fn  void freeMixedLinearComplementarityProblem(MixedLinearComplementarityProblem* problem)
   *  \brief function to delete a MixedLinearComplementarityProblem
   *  \param problem  pointer to a MixedLinearComplementarityProblem to delete
   */
  void freeMixedLinearComplementarityProblem(MixedLinearComplementarityProblem* problem);

  /** create empty MLCP
   * \return empy MLCP
   */
  MixedLinearComplementarityProblem* newMLCP(void);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

