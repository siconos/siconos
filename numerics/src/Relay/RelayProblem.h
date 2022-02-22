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
#ifndef RELAY_PROBLEM_H
#define RELAY_PROBLEM_H

/*!\file RelayProblem.h
  \brief Structure used to define a Relay (dual or primal) Problem

*/

#include <stdio.h>        // for FILE
#include "NumericsFwd.h"  // for RelayProblem, NumericsMatrix
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

/** \struct RelayProblem RelayProblem.h
 * \brief Struct defining a Relay problem
 */
struct RelayProblem
{
  int size;          /**< size dim of the problem */
  NumericsMatrix* M; /**< M matrix of the Relay */
  double* q;        /**< q vector */
  double* lb;       /**< lb upper bound */
  double* ub;       /**< ub lower bound */
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** displays the problem onto the screen 
   * \param[in] p Relay problem to be displayed
   */
  void Relay_display(RelayProblem* p);
  
  /** Saves relay struct into a file.
      \param[in] problem structure
      \param[in] file, file descriptor
      \return file status (1 if everything has worked properly)
   */
  int relay_printInFile(RelayProblem*  problem, FILE* file);

  /* create a new Relay Problem (ptr) and set structure member to null.
   * \return a pointer to a relay 
   */
  RelayProblem* relayProblem_new(void);
  
  /** read a RelayProblem from a file descriptor
   * \param file descriptor
   * \return problem the problem to read
   */
  RelayProblem*  relay_newFromFile(FILE* file);

  /** read a RelayProblem from a file (.dat or hdf5 if fclib is on) from its filename
   * \param filename the name of the input file
   * \return problem the problem to read
   */
  RelayProblem* relay_new_from_filename(const char * filename);

  /** Release memory for the problem structure
      \param[inout] problem, Relay problem structure to be freed.
  */
  void freeRelay_problem(RelayProblem* problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
