/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef RELAY_PROBLEM_H
#define RELAY_PROBLEM_H


#include <assert.h>
/*! \page RelayProblem Relay or box-constrained AVI problems
  \section relayIntro The problem
  Find \f$(z,w)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  w = M z + q\\
  -z \in \mathcal{N}_{K}(w)\\
  \end{array},
  \right.
  \f$

  where M is an (\f$ n \times n \f$)-matrix, q, z and w are n-dimensional vectors, K is the box
  defined by \f$K=\{x\in\mathbb{R}^n \mid lb_i \leq x_i \leq ub_i, i = 1, ..., n \}\f$ and
  \f$\mathcal{N}\f$ is the normal cone to K.

  The solvers and their parameters are described in \ref RelaySolvers. \n

  \section relaySolversList Available solvers

  - relay_avi_caoferris based on an algorithm by Cao and Ferris
  - relay_path using the PATH solver

  Using an LCP reformulation (splitting z in positive and negative part), we have the following
  available solvers:

  - relay_enum which solves the LCP using the enumerative method
  - relay_lexicolemke which solves the LCP using Lemke's algorithm

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
