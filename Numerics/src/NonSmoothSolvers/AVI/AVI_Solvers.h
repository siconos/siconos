/* Siconos-Numerics, Copyright INRIA 2005-2014.
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
#ifndef AVI_SOLVERS_H
#define AVI_SOLVERS_H

/*!\file AVI_Solvers.h
  \brief Subroutines for the resolution of Affine Variational Inequalities.\n

  \author siconos-team@lists.gforge.inria.fr
*/

/*! \page AVISolvers Affine Variational Inequalities Solvers

This page gives an overview of the available solvers for AVI and their required parameters.

For each solver, the input argument are:
- an AffineVariationalInequalities
- the unknown z
- the value of the function F(z)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions struct, which contains iparam and dparam

\section aviCaoFerris The Cao-Ferris algorithm

Direct solver for AVI based on pivoting method principle for (degenerated) problem.

 function: avi_caoferris() \n
 parameters:
- iparam[0] (in): max. number of iterations
- iparam[1] (in,out): boolean, 1 if memory has been allocated, 0 otherwise
- iparam[1] (out): number of iterations processed
*/

#include "AffineVariationalInequalities.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** set the default solver parameters and perform memory allocation for an AVI resolution
      \param[in] problem the AffineVariationalInequalities struct which handles the AVI
      \param options the pointer to the array of options to set
      \param solverId the solver to invoke
      \return info termination value
  */
  int avi_setDefaultSolverOptions(AffineVariationalInequalities* problem, SolverOptions* options, int solverId);


  /** avi_caoferris is a direct solver for AVI based on pivoting method principle for degenerate problem \n
   * Choice of pivot variable is performed via lexicographic ordering \n
   *  Ref: "A Pivotal Method for Affine Variational Inequalities" Menglin Cao et Michael Ferris (1996)\n
   * \param[in] problem structure that represents the AVI (M, q, K)
   * \param[in,out] z on call contains the initial solution and on return holds the solution of the problem.
   * \param[in,out] w defined as Mz + q
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   * \return info about the convergence: 0 ok; 1 ...
   *\author Olivier Huber
   */
  int avi_caoferris(AffineVariationalInequalities* problem, double *z, double* w, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int avi_caoferris_setDefaultSolverOptions(SolverOptions* options);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution \n
   * of the AVI : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * \param[in] problem structure that represents the AVI (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[in] tolerance threshold used to validate the solution: if the error
   * is less than this value, the solution is accepted
   * \param[in,out] error
   * \return status: 0 : convergence, 1: error > tolerance
   * \author Olivier Huber
   */
  //int avi_compute_error(AffineVariationalInequalities* problem, double *z , double *w, double tolerance, double* error);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution \n
  * of the AVI : \n
  * \f$
  *    0 \le z \perp Mz + q \ge 0
  * \f$
  * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
  * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
  * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
  * \param[in] n size of the AVI
  * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
  * \param[in,out] error
  * \author Olivier Huber
  */
//  void avi_compute_error_only(int n,  double *z , double *w, double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

