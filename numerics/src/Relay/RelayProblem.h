/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef RELAY_PROBLEM_H
#define RELAY_PROBLEM_H


#include <assert.h>
/*! \page RelayProblem Relay or box-constrained AVI problems
  \section relayIntro The problem
  Find \f$(z,w)\f$ such that:
  \f{equation*}
  \left\lbrace
  \begin{array}{l}
  w = M z + q\\
  -w \in \mathcal{N}_{K}(z)\\
  \end{array},
  \right.
  \f}
  where M is an (\f$ n \times n \f$)-matrix, q, z and w are n-dimensional vectors, K is the box
  defined by \f$K=\{x\in\mathbb{R}^n \mid lb_i \leq x_i \leq ub_i, i = 1, ..., n \}\f$ and
  \f$\mathcal{N}_K(z)\f$ is the normal cone to \f$K\f$ at \f$z\f$.

  The solvers and their parameters are described in \ref RelaySolvers. \n

  \section relaySolversList Available solvers

  The "direct" solvers are
  - relay_avi_caoferris() based on an algorithm by Cao and Ferris for AVI with a polytopic set \f$K\f$.
  - relay_path() using the PATH solver

  Using an LCP reformulation (splitting z in positive and negative part), we have the following
  available solvers:

  - relay_enum() which solves the LCP using the enumerative method
  - relay_lexicolemke() which solves the LCP using Lemke's algorithm

  (see the functions/solvers list in Relay_Solvers.h)

*/

/*!\file RelayProblem.h
  \brief Structure used to define a Relay (dual or primal) Problem

  \author Franck Perignon
*/

#include "NumericsMatrix.h"

/** \struct RelayProblem RelayProblem.h
 * \brief Struct defining a Relay problem
 */
typedef struct
{
  int size;          /**< size dim of the problem */
  NumericsMatrix* M; /**< M matrix of the Relay */
  double* q;        /**< q vector */
  double* lb;       /**< lb upper bound */
  double* ub;       /**< ub lower bound */
} RelayProblem;

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
