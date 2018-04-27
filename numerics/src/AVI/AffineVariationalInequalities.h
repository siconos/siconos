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
#ifndef AVI_PROBLEM_H
#define AVI_PROBLEM_H

/*!\file AffineVariationalInequalities.h
 * \brief Definitions for AVI
 *
 * \author Olivier Huber
*/

#include "NumericsFwd.h"
#include <stdio.h>

#include "SiconosSets.h"

/** \struct AffineVariationalInequalities AffineVariationalInequalities.h
 *  \brief Structure that contains and defines an \ref AVI
 *
 *   The problem is the following: given a matrix \f$M\f$ and \f$q\f$, find \f$z\f$ such that:
 *   \f{equation*}{
 *   \langle x - z, q + Mz \rangle \geq 0\quad \text{for all }x\in K.
 *   \f}
 */
struct AffineVariationalInequalities
{
  unsigned size;     /**< size of the problem */
  NumericsMatrix* M; /**< M matrix of the AVI (see the mathematical description)*/
  double* q;         /**< vector of the AVI (see the mathematical description)*/
  double* d;         /**< Covering vector (optional) */
  polyhedron_set poly;  /**< Polyhedra where the solution has to belong */
  double* lb;        /**< Lower bounds for the variables */
  double* ub;        /**< Upper bounds for the variables */
  void* cones;       /**< Non-oyhedral Cones where the variable lives (not implemented yet) */
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** \fn void AVI_display(AffineVariationalInequalities* problem)
   *  \brief function to display a AffineVariationalInequalities
   *  \param  avi pointer to a AffineVariationalInequalities to display
   */
  void AVI_display(AffineVariationalInequalities* avi);

  /** \fn int AVI_printInFile(AffineVariationalInequalities*  problem, FILE* file)
   *  \brief function to write in a file a AffineVariationalInequalities
   *  \param avi pointer to a AffineVariationalInequalities to print
   *  \param file pointer to a FILE
   *  \return 1 if successfull
   */
  int AVI_printInFile(AffineVariationalInequalities* avi, FILE* file);

  /** \fn  int AVI_newFromFile(AffineVariationalInequalities* avi, FILE* file)
   *  \brief function to read and create a AffineVariationalInequalities
   *   from a file
   *  \param avi pointer to a AffineVariationalInequalities to create
   *  \param file pointer to a FILE
   *  \return 1 if successfull
   */
  int AVI_newFromFile(AffineVariationalInequalities* avi, FILE* file);

  /** \fn  int AVI_newFromFilename(AffineVariationalInequalities* avi, FILE* file)
   *  \brief function to read and create a AffineVariationalInequalities
   *   from a file
   *  \param avi pointer to a AffineVariationalInequalities to create
   *  \param filename that contains the AVI
   *  \return 1 if successfull
   */
  int AVI_newFromFilename(AffineVariationalInequalities* avi, char* filename);

  /** \fn  void freeAVI(AffineVariationalInequalities* avi)
   *  \brief function to delete a AffineVariationalInequalities
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

