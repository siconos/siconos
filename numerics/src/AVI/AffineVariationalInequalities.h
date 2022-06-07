/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef AVI_PROBLEM_H
#define AVI_PROBLEM_H

/*!\file AffineVariationalInequalities.h
 * \brief Definitions for AVI
 *
*/

#include <stdio.h>        // for FILE
#include "NumericsFwd.h"  // for AffineVariationalInequalities, NumericsMatrix
#include "SiconosSets.h"  // for polyhedron_set

/** Structure that contains and defines an AVI

  The problem is : given a matrix \f$M\f$ and a vector\f$q\f$, find \f$z\f$ such that

 \rst
 
 .. math::

     \langle x - z, q + Mz \rangle \geq 0 \ \text{for all }x\in K


See  :ref:`avi_problem`.
 
\endrst


 */
struct AffineVariationalInequalities
{
  size_t size;     /**< size of the problem */
  NumericsMatrix* M; /**< M matrix of the AVI (see the mathematical description)*/
  double* q;         /**< vector of the AVI (see the mathematical description)*/
  double* d;         /**< Covering vector (optional) */
  polyhedron_set poly;  /**< Polyhedra where the solution has to belong */
  double* lb;        /**< Lower bounds for the variables */
  double* ub;        /**< Upper bounds for the variables */
  void* cones;       /**< Non-polyhedral Cones where the variable lives (not implemented yet) */
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Affine Variational Inequalities display
   *  \param avi pointer to the AffineVariationalInequalities to display
   */
  void AVI_display(AffineVariationalInequalities* avi);

  /** write AVI to file
   *  \param avi pointer to a AffineVariationalInequalities to print
   *  \param file pointer to a FILE
   *  \return 1 if successfull
   */
  int AVI_printInFile(AffineVariationalInequalities* avi, FILE* file);

  /** read from file and create AffineVariationalInequalities
   *  \param avi pointer to a AffineVariationalInequalities to create
   *  \param file pointer to a FILE
   *  \return 1 if successfull
   */
  int AVI_newFromFile(AffineVariationalInequalities* avi, FILE* file);

  /** function to read and create a AffineVariationalInequalities
   *   from a file
   *  \param avi pointer to a AffineVariationalInequalities to create
   *  \param filename that contains the AVI
   *  \return 1 if successfull
   */
  int AVI_newFromFilename(AffineVariationalInequalities* avi, char* filename);

  /** function to delete a AffineVariationalInequalities
   *  \param avi  pointer to a AffineVariationalInequalities to delete
   */
  void freeAVI(AffineVariationalInequalities* avi);

  /** Create an empty AVI struct
   * \return an empty AffineVariationalInequalities*/
  AffineVariationalInequalities* newAVI(void);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

