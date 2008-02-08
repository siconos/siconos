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
#ifndef LCP_SOLVERS_H
#define LCP_SOLVERS_H

/*!\file LCP_Solvers.h
  \brief Subroutines for the resolution of Linear Complementarity Problems.\n

  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

#include "Numerics_Options.h"
#include "LinearComplementarity_Problem.h"
#include "Solver_Options.h"

#ifdef __cplusplus
extern "C" {
#endif

  /** lcp_qp uses a quadratic programm formulation for solving a LCP
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   *                0 : convergence  / minimization sucessfull\n
   *                1 : Too Many iterations\n
   *              2 : Accuracy insuficient to satisfy convergence criterion\n
   *                5 : Length of working array insufficient\n
   *                Other : The constraints are inconstent\n
   * \param[in-out] options, structure used to define the solver and its parameters.
   *
   * \author Vincent Acary
   */
  void lcp_qp(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_cpg is a CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   0: convergence\n
   1: iter = itermax\n
   2: negative diagonal term\n
   3: pWp nul
   * \param[in-out] options, structure used to define the solver and its parameters.
   *

   \author Mathieu Renouf.
  */
  void lcp_cpg(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for LCP.\n
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options, structure used to define the solver and its parameters.
   \author Mathieu Renouf
  */
  void lcp_pgs(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_rpgs (Regularized Projected Gauss-Seidel ) is a solver for LCP, able to handle matrices with null diagonal terms.\n
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   * \param[in-out] options, structure used to define the solver and its parameters.
   *

   \author Mathieu Renouf & Pascal Denoyelle
   \todo Sizing the regularization paramter and apply it only on null diagnal term

  */
  void lcp_rpgs(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_psor Projected Succesive over relaxation solver for LCP. See cottle, Pang Stone Chap 5
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   * \param[in-out] options, structure used to define the solver and its parameters.

   \author  Vincent Acary

   \todo use the relax parameter
   \todo add test
   \todo add a vector of relaxation parameter wtith an auto sizing (see SOR algorithm for linear solver.)

  */
  void lcp_psor(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_nsqp use a quadratic programm formulation for solving an non symmetric LCP
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   *                0 : convergence  / minimization sucessfull\n
   *                1 : Too Many iterations\n
   *            2 : Accuracy insuficient to satisfy convergence criterion\n
   *                5 : Length of working array insufficient\n
   *                Other : The constraints are inconstent\n
   * \param[in-out] options, structure used to define the solver and its parameters.
   *
   * \author Vincent Acary
   *
   */
  void lcp_nsqp(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_latin (LArge Time INcrements) is a basic latin solver for LCP.
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : Cholesky Factorization failed \n
   3 : nul diagonal term\n
   * \param[in-out] options, structure used to define the solver and its parameters.

   \author Nineb Sheherazade.
  */
  void lcp_latin(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_latin_w (LArge Time INcrements) is a basic latin solver with relaxation for LCP.\n
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : Cholesky Factorization failed \n
   3 : nul diagonal term\n
   * \param[in-out] options, structure used to define the solver and its parameters.
   *

   \author Nineb Sheherazade.
  */
  void lcp_latin_w(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_lexicolemke is a direct solver for LCP based on pivoting method principle for degenerate problem \n
      Choice of pivot variable is performed via lexicographic ordering \n
      Ref: "The Linear Complementarity Problem" Cottle, Pang, Stone (1992)\n
      * \param[in] problem, structure that represents the LCP (M, q...)
      * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
      * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
      * \param[out] info, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : negative diagonal term\n
      * \param[in-out] options, structure used to define the solver and its parameters.
      *
      \param dparamLCP  On enter/return, a vetor of doubles (not used).\n


      \author Mathieu Renouf

  */
  void lcp_lexicolemke(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_newton_min uses a nonsmooth Newton method based on the min formulation  (or max formulation) of the LCP
      \f$
      0 \le z \perp w \ge 0 \Longrightarrow \min(w,\rho z)=0 \Longrightarrow w = \max(0,w - \rho z)
      \f$

      \f$
      H(z) = H(\left[ \begin{array}{c} z \\ w \end{array}\right])= \left[ \begin{array}{c} w-Mz-q \\ min(w,\rho z) \end{array}\right] =0\\
      \f$


      References: Alart & Curnier 1990, Pang 1990
      * \param[in] problem, structure that represents the LCP (M, q...)
      * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
      * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
      * \param[out] info, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : Problem in resolution in DGESV\n
      *                0 : convergence  / minimization sucessfull\n
      *                1 : Too Many iterations\n
      *               2 : Accuracy insuficient to satisfy convergence criterion\n
      *                5 : Length of working array insufficient\n
      *                Other : The constraints are inconstent\n
      * \param[in-out] options, structure used to define the solver and its parameters.
      *
      \author Vincent Acary

      \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
      \todo Add rules for the computation of the penalization rho
      \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
  */
  void lcp_newton_min(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** lcp_newton_FB use a nonsmooth newton method based on the Fischer-Bursmeister convex function
   *
   *
   * \f$
   *   0 \le z \perp w \ge 0 \Longrightarrow \phi(z,w)=\sqrt{z^2+w^2}-(z+w)=0
   * \f$

   * \f$
   *   \Phi(z) = \left[ \begin{array}{c}  \phi(z_1,w_1) \\ \phi(z_1,w_1) \\ \vdots \\  \phi(z_n,w_n)  \end{array}\right] =0\\
   * \f$
   *
   *
   * References: Alart & Curnier 1990, Pang 1990
   *
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *                2 - negative diagonal term\n
   *
   * \param[in-out] options, structure used to define the solver and its parameters.
   * \author Vincent Acary
   *
   * \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
   * \todo Add rules for the computation of the penalization rho
   * \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
   */
  void lcp_newton_FB(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /**
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in-out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[out] info, an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   * \param[in-out] options, structure used to define the solver and its parameters.

   \author Olivier Bonnefon
  */
  void lcp_path(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** This function checks the validity of the vector z as a solution \n
   * of the LCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
   * \param[in] problem, structure that represents the LCP (M, q...)
   * \param[in,out] z, a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w, a n-vector of doubles which returns the solution of the problem.
   * \param[in] tolerance
   * \return status: 0 : convergence, 1: error > tolerance
   * \author Pascal Denoyelle
   */
  int filter_result_LCP(LinearComplementarity_Problem* problem, double *z , double *w, double tolerance);

  /**
   * This function checks the validity of the vector z as a solution \n
   * of the LCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
   * \author Vincent Acary form the routine  filter_result_LCP.c of Pascal Denoyelle
   */
  int lcp_compute_error(int n, double *M , double *q , double *z , int chat, double *w, double * error);

  /** This function checks the validity of the vector z as a solution \n
   * of the LCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * where M is a sparse block matrix defined in argument blmat.\n
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
   * \author Pascal Denoyelle
   */
  int filter_result_LCP_block(SparseBlockStructuredMatrix *blmat, double *q , double *z , double tol, int chat, double *w);

  /*   /\** Function used to extract from LCP matrix the part which corresponds to non null z */
  /*    *\/ */
  /*   int extractLCP( NumericsMatrix* MGlobal, double *z , int *indic, int *indicop, double *submatlcp , double *submatlcpop, */
  /*      int *ipiv , int *sizesublcp , int *sizesublcpop); */

  /*   /\** Function used to solve the problem with sub-matrices from extractLCP */
  /*    *\/ */
  /*   int predictLCP(int sizeLCP, double* q, double *z , double *w , double tol, */
  /*    int *indic , int *indicop , double *submatlcp , double *submatlcpop , */
  /*     int *ipiv , int *sizesublcp , int *sizesublcpop , double *subq , double *bufz , double *newz); */

#ifdef __cplusplus
}
#endif

#endif

