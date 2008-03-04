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
#ifndef MLCP_SOLVERS_H
#define MLCP_SOLVERS_H

/*!\file MLCP_Solvers.h
  Solvers for Mixed Linear Complementary Problems (MCLP)
  \author Vincent Acary
*/

/*! \page MLCPSolvers Mixed Linear Complementary Problems Solvers

This page gives an overview of the available solvers for MLCP and their required parameters.

For each solver, the input argument are:
- a MixedLinearComplementarity_Problem
- the unknowns (z,w)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a Solver_Options structure, which handles iparam and dparam

\section mlcpPGS PGS Solver
Projected Gauss-Seidel solver

\bf function: mlcp_pgs() \n
\bf parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error

\section mlcpRPGS RPGS Solver
Regularized Projected Gauss-Seidel, solver for MLCP, able to handle with matrices with null diagonal terms

\bf function: mlcp_rpgs() \n
\bf parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): rho

\section mlcpPSOR PSOR Solver
Projected Succesive over relaxation solver for MLCP. See cottle, Pang Stone Chap 5
\bf function: mlcp_psor() \n
\bf parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): omega

\section mlcpRPSOR RPSOR Solver
Regularized Projected Succesive over relaxation solver for MLCP
\bf function: mlcp_rpsor() \n
\bf parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): omega
- dparam[3] (in): rho

\section mlcpPath Path (Ferris) Solver

\bf function: mlcp_path() \n
\bf parameters:
- dparam[0] (in): tolerance

*/

#include "Numerics_Options.h"
#include "Solver_Options.h"
#include "MixedLinearComplementarity_Problem.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**  mlcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_pgs(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /**  mlcp_rpgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_rpgs(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** mlcp_psor (projected successive overrelaxation method) is a solver for MLCP.
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_psor(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** mlcp_rpsor (regularized projected successive overrelaxation method) is a solver for MLCP.
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_rpsor(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** path solver
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_path(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /**
   * This function checks the validity of the vector z as a solution \n
   * of the MLCP : \n
   * \f$
   *  \left\lbrace
   *   \begin{array}{l}
   *   A u + Cv +a =0\\
   *   D u + Bv +b = w
   *   0 \le v \perp  w \ge 0\\
   *   \end{array}
   *  \right.
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
   * \param[in] problem structure that represents the MLCP (M, q... or (A,B,C...))
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[in] tolerance
   * \param[in,out] error
   * \return status: 0 : convergence, 1: error > tolerance
   * \author Vincent Acary form the routine  filter_result_LCP.c of Pascal Denoyelle
   */
  int mlcp_compute_error(MixedLinearComplementarity_Problem* problem, double *z, double *w, double tolerance, double * error);

#ifdef __cplusplus
}
#endif

#endif
