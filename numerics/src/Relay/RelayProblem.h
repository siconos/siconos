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
#ifndef RELAY_PROBLEM_H
#define RELAY_PROBLEM_H


#include <assert.h>
#include <stdio.h>

/*!\file RelayProblem.h
  \brief Structure used to define a Relay (dual or primal) Problem

  \author Franck Perignon
*/

#include "NumericsFwd.h"
#include "SiconosConfig.h"

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

  /** Relay_display displays on screen a Relay_problem
  * \param[in] p Relay_problem to be displayed
  * \author Vincent Acary
  */
  void Relay_display(RelayProblem* p);

  int relay_printInFile(RelayProblem*  problem, FILE* file);

  int relay_newFromFile(RelayProblem* problem, FILE* file);

  void freeRelay_problem(RelayProblem* problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
