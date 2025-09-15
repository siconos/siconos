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
#ifndef MLCP_PROBLEM_H
#define MLCP_PROBLEM_H

/*!\file MixedLinearComplementarityProblem.h
  \brief Structure used to define a Mixed Linear Complementarity Problem

*/

#include <stdio.h>  // for FILE

#include "NumericsFwd.h"    // for MixedLinearComplementarityProblem, Numerics...
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

/** The Structure that contains and defines MLCProblem. Find \f$ (z,w) \f$ such
 *  that:
 *
 *  \f[ \left\{ \begin{array}{l}
 *   M \ z + q = w \\
 *  w_1=0 \\
 *  0 \le w_{2} \perp v \ge 0
 *  \end{array}
 *  \right.
 *  \text{ with }
 *  z=
 *  \left[
 *  \begin{array}{c}
 *  u\\
 *  v\\
 *  \end{array}
 *  \right]
 *  \text{ and }
 *  w=
 *  \left[
 *  \begin{array}{c}
 *  w_{1}\\
 *  w_{2}\\
 *  \end{array}
 *  \right]
 *  \f]
 *
 *  \f$ u, w_{1} \f$ are vectors of size n.
 *  \f$ v, w_{2} \f$ are vectors of size m.
 *
 *  ABCD format (see  "R. W. {Cottle} and J. {Pang} and R. E. {Stone}", "The
 *  Linear Complementarity Problem, Academic Press, Inc., 1992, Section 1.5 )
 *
 *  \f[
 *  \left[
 *  \begin{array}{cc}
 *  A & C \\
 *  D & B \\
 *  \end{array}
 *  \right]
 *  \f]
 */
struct MixedLinearComplementarityProblem {
  int isStorageType1; /**< boolean for storageType1 1 if the problem
                         is saved using (M,q),  0 otherwise */
  int isStorageType2; /**< boolean for storageType2 1 if the problem
                         is saved using (A,B,C,D,a,b), 0 otherwise*/
  int n;              /**< number of equality constraints */
  int m;              /**< number of complementarity constraints */
  int *blocksRows;    /**< The rows from blocksRows[i] to blocksRows[i+1]-1
                         forms a block of equalities iif bloksIsComp[i]=0,
                         else the block is a complementarity block.
                         The number of total blocks is given by NbBlocks
                         such that blocksRows[NbBlocks] = n+m */
  int *blocksIsComp;  /**< if bloksIsComp[i]=0, then block i formed by the rows
                         from blocksRows[i] to blocksRows[i+1]-1 is an equality
                         block  else the block is a complementarity block.
                      */
  NumericsMatrix *M;  /**< M matrix of the MLCP */
  double *q;          /**< q vector of the MLCP */
  // NumericsMatrix* Bblock; Bblock  ?*/
  double *A; /**< A matrix of the MLCP */
  double *B; /**< B matrix of the MLCP */
  double *C; /**< C matrix of the MLCP */
  double *D; /**< D matrix of the MLCP */
  double *a; /**< a vector of the MLCP */
  double *b; /**< b vector of the MLCP */
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/** function to delete a MixedLinearComplementarityProblem
 *
 *  \param problem  pointer to a MixedLinearComplementarityProblem to delete
 */
void mixedLinearComplementarity_free(MixedLinearComplementarityProblem *problem);

/** create empty MLCP
 *
 *  \return empy MLCP
 */
MixedLinearComplementarityProblem *mixedLinearComplementarity_new(void);

/** display a MLCP
 */
void mixedLinearComplementarity_display(MixedLinearComplementarityProblem *p);

/** function to write in a file a MixedLinearComplementarityProblem
 *
 *  \param problem pointer to a MixedLinearComplementarityProblem to print
 *  \param file pointer to a FILE
 *  \return 0 if ok
 */
int mixedLinearComplementarity_printInFile(MixedLinearComplementarityProblem *problem,
                                           FILE *file);

/** Function to read and create a MixedLinearComplementarityProblem
 *  from a file
 *
 *  \param problem pointer to a MixedLinearComplementarityProblem to create
 *  \param file pointer to a FILE
 *  \return 0 if ok
 */
int mixedLinearComplementarity_newFromFile(MixedLinearComplementarityProblem *problem,
                                           FILE *file);

/** Function to read and create a MixedLinearComplementarityProblem
 *  from a file
 *
 *  \param problem pointer to a MixedLinearComplementarityProblem to create
 *  \param file pointer to a FILE
 *  \return 0 if ok
 */
int mixedLinearComplementarity_newFromFileOld(MixedLinearComplementarityProblem *problem,
                                              FILE *file);

/** Function to read and create a MixedLinearComplementarityProblem
 *  from a file
 *
 *  \param problem pointer to a MixedLinearComplementarityProblem to create
 *  \param filename that contains the mlcp
 *  \return 0 if ok
 */
int mixedLinearComplementarity_newFromFilename(MixedLinearComplementarityProblem *problem,
                                               const char *filename);

/** Function to create a MLCP with ABCD format from M formatted MLCP
 */

MixedLinearComplementarityProblem *mixedLinearComplementarity_fromMtoABCD(
    MixedLinearComplementarityProblem *problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
