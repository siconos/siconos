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
#ifndef LCP_PROBLEM_H
#define LCP_PROBLEM_H

/*!\file LinearComplementarityProblem.h
*/

/*! \page LCProblem Linear Complementarity problems (LCP)
 * \section lcpIntro The problem
 *  The Linear Complementarity problems (LCP) is defined by
 *
 *   Find \f$(z,w)\f$ such that:\n
 *   \f{equation*}{
 *   \begin{cases}
 *   M \ z + q = w \\
 *   0 \le w \perp z \ge 0 \\
 *   \end{cases}
 *   \f}
 *
 * where \f$ w, z, q\f$ are vectors of size \f$n\f$ and \f$ M \f$ is a \f$n\times n\f$ matrix.
 *
 * The notation \f$x \perp y\f$ means that \f$x^Ty =0\f$. Inequalities involving vectors
 * are understood to hold component-wise.
 *
 * From more details on theory and analysis of LCP, we refer to
 *
 * R.W. Cottle, J.S. Pang, and R.E. Stone. <i>The Linear Complementarity Problem.</i> Academic Press, Inc., Boston, MA, 1992.
 *
 *  \section lcpSolversList Available solvers

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
  - lcp_nsgs_SBM(), Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.

  (see also the functions/solvers list in LCP_Solvers.h and numbering in lcp_cst.h)

*/

#include "NumericsMatrix.h"

/** \struct LinearComplementarityProblem LinearComplementarityProblem.h
 *  \brief Structure that contains and defines  \ref LCProblem
 *
 *   Find \f$(z,w)\f$ such that:\n
 *   \f{equation*}{
 *   \begin{cases}
 *   M \ z + q = w \\
 *   0 \le w \perp z \ge 0 \\
 *   \end{cases}
 *   \f}
 *
 * where \f$ w, z, q\f$ are vectors of size \f$n\f$ and \f$ M \f$ is a \f$n\times n\f$ matrix.
 * See \ref LCProblem for more details.
 */
typedef struct
{

  int size; /**<  size of the problem */
  NumericsMatrix* M ;/**< M matrix of the LCP (see the mathematical description)*/
  double * q;/**< vector of the LCP (see the mathematical description)*/
} LinearComplementarityProblem;

#ifdef __cplusplus
extern "C"
{
#endif
  /** \fn void LinearComplementarity_display(LinearComplementarityProblem* problem)
   *  \brief function to display a LinearComplementarityProblem
   *  \param  problem pointer to a LinearComplementarityProblem to display
   */
  void LinearComplementarity_display(LinearComplementarityProblem* problem);

  /** \fn int linearComplementarity_printInFile(LinearComplementarityProblem*  problem, FILE* file)
   *  \brief function to write in a file a LinearComplementarityProblem
   *  \param problem pointer to a LinearComplementarityProblem to print
   *  \param file pointer to a FILE
   */
  int linearComplementarity_printInFile(LinearComplementarityProblem*  problem, FILE* file);

  /** \fn  int linearComplementarity_newFromFile(LinearComplementarityProblem* problem, FILE* file)
   *  \brief function to read and create a LinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a LinearComplementarityProblem to create
   *  \param file pointer to a FILE
   */
  int linearComplementarity_newFromFile(LinearComplementarityProblem* problem, FILE* file);

  /** \fn  void freeLinearComplementarityProblem(LinearComplementarityProblem* problem)
   *  \brief function to delete a LinearComplementarityProblem
   *  \param problem  pointer to a LinearComplementarityProblem to delete
   */
  void freeLinearComplementarityProblem(LinearComplementarityProblem* problem);
#ifdef __cplusplus
}
#endif

#endif

