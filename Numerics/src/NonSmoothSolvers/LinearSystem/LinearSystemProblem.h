/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#ifndef LINEARSYSTEM_PROBLEM_H
#define LINEARSYSTEM_PROBLEM_H

/*!\file LinearComplementarityProblem.h
  \brief Structure used to define a Linear Complementarity Problem

  \author Franck Perignon
*/

/*! \page LCProblem System problems
  \section lcpIntro The problem
  Find \f$(z,w)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ z + q = w =0\\
  \end{array}
  \right.
  \f$

  \f$ w, z, q\f$ are vectors of size n and \f$ M \f$ is a nXn matrix.

  \section lcpSolversList Available solvers

  Use the generic function LinearSystem_driver() to call one the the specific solvers listed below:



*/

#include "NumericsMatrix.h"

/** Linear Complementarity Problem elements
    \param size dim of the problem
    \param M matrix of the linear system
    \param q vector
 */
typedef struct
{
  int size;
  NumericsMatrix* M;
  double * q;
} LinearSystemProblem;

#ifdef __cplusplus
extern "C" {
#endif
  void displayLS(LinearSystemProblem* p);
#ifdef __cplusplus
}
#endif

#endif

