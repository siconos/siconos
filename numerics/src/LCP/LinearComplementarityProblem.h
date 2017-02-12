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
#ifndef LCP_PROBLEM_H
#define LCP_PROBLEM_H

/*!\file LinearComplementarityProblem.h
*/

/*! \page LCProblem Linear Complementarity problems (LCP)
 * \section lcpIntro The problem
 *  The Linear Complementarity problem (LCP) is defined by
 *
 *   Find \f$(z,w)\f$ such that:
 *   \f{equation*}{
 *   \begin{cases}
 *   M \ z + q = w \\
 *   0 \le w \perp z \ge 0
 *   \end{cases},
 *   \f}
 * where \f$ w, z, q\f$ are vectors of size \f$n\f$ and \f$ M \f$ is a \f$n\times n\f$ matrix.
 *
 * The notation \f$x \perp y\f$ means that \f$x^Ty =0\f$. Inequalities involving vectors
 * are understood to hold component-wise.
 *
 * From more details on theory and analysis of LCP, we refer to
 *
 * R.W. Cottle, J.S. Pang, and R.E. Stone. <i>The Linear Complementarity Problem.</i> Academic Press, Inc., Boston, MA, 1992.
 *
* The problem is stored and given to the solver in numerics thanks to
   the C structure LinearComplementarityProblem.
 
 *  \section lcpSolversList Available solvers

  Use the generic functions lcp_driver_DenseMatrix() to call one the the specific solvers listed below:

  \subsection lcpDirectSolvers Direct solvers
  - lcp_lexicolemke(), direct solver for LCP based on pivoting method principle for degenerate problem:
  the choice of pivot variable is performed via lexicographic ordering.
  - lcp_pivot(), generic solver for pivot-based methods: Bard, Murty and Lemke rules are implemented.
  - lcp_enum(), enumerative solver (brute-force method which tries every possible solution).
  - lcp_path(), interface to the PATH solver

  \subsection lcpIterativeSolvers Iterative solvers
  - lcp_cpg(), CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.
  - lcp_pgs(), PGS is a basic Projected Gauss-Seidel solver for LCP.
  - lcp_rpgs(), Regularized Projected Gauss-Seidel, is a solver for LCP, able to handle matrices with null diagonal terms.
  - lcp_psor(), Projected Succesive over relaxation solver for LCP. See Cottle, Pang Stone Chap 5.
  - lcp_latin(), (LArge Time INcrements) is a basic latin solver for LCP.
  - lcp_latin_w(), (LArge Time INcrements) is a basic latin solver with relaxation for LCP.
  - lcp_nsgs_SBM(), Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.

  \subsection lcpEquationBasedSolvers Equation-based solvers
  - lcp_newton_min(), nonsmooth Newton method based on the min formulation of the LCP.
  - lcp_newton_FB(), uses a nonsmooth newton method based on the Fischer-Bursmeister NCP function.
  - lcp_newton_minFB(), nonsmooth Newton method combining the min and FB functions.

  \subsection lcpReformulation QP-reformulation
  - lcp_qp(), quadratic programm formulation
  - lcp_nsqp(), quadratic programm formulation for solving an non symmetric LCP

  (see also the functions/solvers list in LCP_Solvers.h and numbering in lcp_cst.h)

*/

#include "NumericsFwd.h"
#include "SiconosConfig.h"

#include <stdio.h>
/** \struct LinearComplementarityProblem LinearComplementarityProblem.h
 *  \brief Structure that contains and defines \ref LCProblem
  */
struct LinearComplementarityProblem
{

  int size; /**<  size of the problem */
  NumericsMatrix* M; /**< M matrix of the LCP (see the mathematical description)*/
  double * q; /**< vector of the LCP (see the mathematical description)*/
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** \fn void linearComplementarity_display(LinearComplementarityProblem* problem)
   *  \brief function to display a LinearComplementarityProblem
   *  \param  problem pointer to a LinearComplementarityProblem to display
   */
  void linearComplementarity_display(LinearComplementarityProblem* problem);

  /** \fn int linearComplementarity_printInFile(LinearComplementarityProblem*  problem, FILE* file)
   *  \brief function to write in a file a LinearComplementarityProblem
   *  \param problem pointer to a LinearComplementarityProblem to print
   *  \param file pointer to a FILE
   *  \return 0 if ok
   */
  int linearComplementarity_printInFile(LinearComplementarityProblem*  problem, FILE* file);

  /** \fn  int linearComplementarity_newFromFile(LinearComplementarityProblem* problem, FILE* file)
   *  \brief function to read and create a LinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a LinearComplementarityProblem to create
   *  \param file pointer to a FILE
   *  \return 0 if ok
   */
  int linearComplementarity_newFromFile(LinearComplementarityProblem* problem, FILE* file);

  /** \fn  int linearComplementarity_newFromFilename(LinearComplementarityProblem* problem, FILE* file)
   *  \brief function to read and create a LinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a LinearComplementarityProblem to create
   *  \param filename that contains the lcp
   *  \return 0 if ok
  */
  int linearComplementarity_newFromFilename(LinearComplementarityProblem* problem, char* filename);

  /** \fn  void freeLinearComplementarityProblem(LinearComplementarityProblem* problem)
   *  \brief function to delete a LinearComplementarityProblem
   *  \param problem  pointer to a LinearComplementarityProblem to delete
   */
  void freeLinearComplementarityProblem(LinearComplementarityProblem* problem);

  /** Create new LCP and clear its fields
   * \return a LinearComplementarityProblem
   */
  LinearComplementarityProblem* newLCP(void);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

