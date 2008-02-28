/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#ifndef MLCP_PROBLEM_H
#define MLCP_PROBLEM_H

/*! \page MLCProblem Mixed Linear Complementarity problems (MLCP)
  \section lcpIntro The problem
  Find \f$(z,w)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ z + q = w \\
  \end{array}
  \right.
  \f$

  \f$
  z=
  \left\lbrace
  \begin{array}{c}
  u\\
  v\\
  \end{array}
  \right.
  \f$

  \f$
  w=
  \left\lbrace
  \begin{array}{c}
  w_{1}\\
  w_{2}\\
  \end{array}
  \right.
  \f$

  \f$ u, w_{1}\f$ are vectors of size n.
  \f$ v, w_{2}\f$ are vectors of size m.

  \f$ w_1=0\f$.
  \f$ 0 \le w_{2} \perp v \ge 0 \f$

  \section mlcpSolversList Available solvers

  The solvers and their parameters are described in \ref MLCPSolvers . \n
  Use the generic function mlcp_driver() to call one the the specific solvers listed below:


  (see the functions/solvers list in MLCP_Solvers.h)

*/

/*!\file MixedLinearComplementarity_Problem.h
  \brief Structure used to define a Mixed Linear Complementarity Problem

*/

#include "NumericsMatrix.h"

/** Mixed Linear Complementarity Problem elements
    \param n dim of the linear constarints
    \param m dim of the Complementarity constraints
    \param M matrix of the MLCP
    \param q vector
 */
typedef struct
{
  int n;
  int m;
  NumericsMatrix* M;
  double * q;
} MixedLinearComplementarity_Problem;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif

