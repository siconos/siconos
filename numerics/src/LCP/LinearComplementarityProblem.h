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
#ifndef LCP_PROBLEM_H
#define LCP_PROBLEM_H

/*!\file LinearComplementarityProblem.h
 */

#include <stdio.h>  // for FILE

#include "NumericsFwd.h"    // for LinearComplementarityProblem, NumericsMatrix
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

/**
   Structure that contains and defines a LCP Problem.
*/
struct LinearComplementarityProblem {
  int size;          /**<  size of the problem */
  NumericsMatrix *M; /**< M matrix of the LCP (see the mathematical description)*/
  double *q;         /**< vector of the LCP (see the mathematical description)*/
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/** display a LinearComplementarityProblem
 *
 *  \param  problem pointer to a LinearComplementarityProblem to display
 */
void linearComplementarity_display(LinearComplementarityProblem *problem);

/** function to write in a file a LinearComplementarityProblem
 *
 *  \param problem pointer to a LinearComplementarityProblem to print
 *  \param file pointer to a FILE
 *  \return 0 if ok
 */
int linearComplementarity_printInFile(LinearComplementarityProblem *problem, FILE *file);

/** Read and create a LinearComplementarityProblem
 *  from a file
 *
 *  \param problem pointer to a LinearComplementarityProblem to create
 *  \param file pointer to a FILE
 *  \return 0 if ok
 */
int linearComplementarity_newFromFile(LinearComplementarityProblem *problem, FILE *file);

/** Read and create a LinearComplementarityProblem
 *  from a file
 *
 *  \param problem pointer to a LinearComplementarityProblem to create
 *  \param filename that contains the lcp
 *  \return 0 if ok
 */
int linearComplementarity_newFromFilename(LinearComplementarityProblem *problem,
                                          const char *filename);

/** Delete a LinearComplementarityProblem
 *
 *  \param problem  pointer to a LinearComplementarityProblem to delete
 */
void freeLinearComplementarityProblem(LinearComplementarityProblem *problem);

/** Create new LCP and clear its fields
 *
 *  \return a LinearComplementarityProblem
 */
LinearComplementarityProblem *newLCP(void);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
