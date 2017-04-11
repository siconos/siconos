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

  /** avi_pathavi is using PATHVI, a direct solver for VI based on pivoting method principle for degenerate problem \n
   *  Ref: "A structure-preserving Pivotal Method for Affine Variational Inequalities" Y. Kim, O. Huber, M.C. Ferris, Math Prog B (2017)\n
   * \param[in] problem structure that represents the AVI (M, q, K)
   * \param[in,out] z on call contains the initial solution and on return holds the solution of the problem.
   * \param[in,out] w defined as Mz + q
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   * \return info about the convergence: 0 ok; 1 ...
   *\author Olivier Huber
   */
  int avi_pathavi(AffineVariationalInequalities* problem, double *z, double *w, SolverOptions* options);

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

