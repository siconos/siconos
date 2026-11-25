/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef LCP_PROBLEM_H
#define LCP_PROBLEM_H

/*!\file LinearComplementarityProblem.h
 */

#include <stdio.h>  // for FILE

#include "NumericsFwd.h"  // for NumericsMatrix

/**
   Structure that contains and defines a LCP Problem.
*/
struct LinearComplementarityProblem {
  int size;          /**< size of the problem */
  NumericsMatrix* M; /**< M matrix of the LCP (see the mathematical description)*/
  double* q;         /**< vector of the LCP (see the mathematical description)*/
};

#if defined(__cplusplus)
extern "C" {
#endif

/** display a LinearComplementarityProblem
 *
 *  \param  problem pointer to a LinearComplementarityProblem to display
 */
void linearComplementarity_display(LinearComplementarityProblem* problem);

/** function to write in a file a LinearComplementarityProblem
 *
 *  \param problem pointer to a LinearComplementarityProblem to print
 *  \param file pointer to a FILE
 *  \return 0 if ok
 */
int linearComplementarity_printInFile(LinearComplementarityProblem* problem, FILE* file);

/** Read and create a LinearComplementarityProblem
 *  from a file
 *
 *  \param problem pointer to a LinearComplementarityProblem to create
 *  \param file pointer to a FILE
 *  \return 0 if ok
 */
int linearComplementarity_newFromFile(LinearComplementarityProblem* problem, FILE* file);

/** Read and create a LinearComplementarityProblem
 *  from a file
 *
 *  \param problem pointer to a LinearComplementarityProblem to create
 *  \param filename that contains the lcp
 *  \return 0 if ok
 */
int linearComplementarity_newFromFilename(LinearComplementarityProblem* problem,
                                          const char* filename);

/** Delete a LinearComplementarityProblem
 *
 *  \param problem  pointer to a LinearComplementarityProblem to delete
 */
void freeLinearComplementarityProblem(LinearComplementarityProblem* problem);

/** Create new LCP and clear its fields
 *
 *  \return a LinearComplementarityProblem
 */
LinearComplementarityProblem* newLCP(void);

#if defined(__cplusplus)
}
#endif

#endif
