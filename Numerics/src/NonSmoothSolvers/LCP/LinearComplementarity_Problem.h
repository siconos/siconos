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
#ifndef LCP_PROBLEM_H
#define LCP_PROBLEM_H

/*!\file LinearComplementarity_Problem.h
  \brief Structure used to define a Linear Complementarity Problem

  \author Franck Perignon
*/

/*! \page LCProblem Linear Complementarity problems (LCP)
  \section lcpIntro The problem
  Find \f$(z,w)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ z + q = w \\
  0 \le w \perp z \ge 0 \\
  \end{array}
  \right.
  \f$

  \f$ w, z, q\f$ are vectors of size n and \f$ M \f$ is a nXn matrix.

  \section lcpSolversList Available solvers

  The solvers and their parameters are described in \ref LCPSolvers . \n
  Use the generic function lcp_driver() to call one the the specific solvers listed below:

  - lcp_qp(), quadratic programm formulation
  - lcp_cpg(), CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.
  - lcp_pgs(), PGS is a basic Projected Gauss-Seidel solver for LCP.
  - lcp_rpgs(), Regularized Projected Gauss-Seidel, is a solver for LCP, able to handle matrices with null diagonal terms
  - lcp_psor(), Projected Succesive over relaxation solver for LCP. See cottle, Pang Stone Chap 5
  - lcp_nsqp(), quadratic programm formulation for solving an non symmetric LCP
  - lcp_latin(), (LArge Time INcrements) is a basic latin solver for LCP.
  - lcp_latin_w(), (LArge Time INcrements) is a basic latin solver with relaxation for LCP
  - lcp_lexicolemke(), direct solver for LCP based on pivoting method principle for degenerate problem.\n
  Choice of pivot variable is performed via lexicographic ordering
  Ref: "The Linear Complementarity Problem" Cottle, Pang, Stone (1992)\n
  - lcp_newton_min(), nonsmooth Newton method based on the min formulation  (or max formulation) of the LCP
  - lcp_newton_FB(), uses a nonsmooth newton method based on the Fischer-Bursmeister convex function
  - lcp_GaussSeidel_SBM(), Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.

  (see the functions/solvers list in LCP_Solvers.h)

*/

#include "NumericsMatrix.h"

/** Linear Complementarity Problem elements
    \param size dim of the problem
    \param M matrix of the LCP
    \param q vector
    \param isComplete equal to 0 if some information is missing or wrong for the problem (M or q = NULL, inconsistent sizes), else equal to 1.
 */
typedef struct
{
  int size;
  NumericsMatrix* M;
  double * q;
} LinearComplementarity_Problem;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif

