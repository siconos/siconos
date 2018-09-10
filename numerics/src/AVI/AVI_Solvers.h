/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
  \brief Subroutines for the resolution of Affine Variational Inequalities.

*/

#include "AffineVariationalInequalities.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** avi_caoferris is a direct solver for AVI based on pivoting method principle for degenerate problem 
   * Choice of pivot variable is performed via lexicographic ordering 
   *  Ref: "A Pivotal Method for Affine Variational Inequalities" Menglin Cao et Michael Ferris (1996)
   * \param[in] problem structure that represents the AVI (M, q, K)
   * \param[in,out] z on call contains the initial solution and on return holds the solution of the problem.
   * \param[in,out] w defined as Mz + q
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   * \return info about the convergence: 0 ok; 1 ...
   */
  int avi_caoferris(AffineVariationalInequalities* problem, double *z, double* w, SolverOptions* options);

  /** avi_pathavi is using PATHVI, a direct solver for VI based on pivoting method principle for degenerate problem 
   *  Ref: "A structure-preserving Pivotal Method for Affine Variational Inequalities" Y. Kim, O. Huber, M.C. Ferris, Math Prog B (2017)
   * \param[in] problem structure that represents the AVI (M, q, K)
   * \param[in,out] z on call contains the initial solution and on return holds the solution of the problem.
   * \param[in,out] w defined as Mz + q
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   * \return info about the convergence: 0 ok; 1 ...
   */
  int avi_pathavi(AffineVariationalInequalities* problem, double *z, double *w, SolverOptions* options);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution 
   * of the AVI : 
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ 
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. 
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.
   * \param[in] problem structure that represents the AVI (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[in] tolerance threshold used to validate the solution: if the error
   * is less than this value, the solution is accepted
   * \param[in,out] error
   * \return status: 0 : convergence, 1: error > tolerance
   */
  //int avi_compute_error(AffineVariationalInequalities* problem, double *z , double *w, double tolerance, double* error);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution 
  * of the AVI : 
  * \f$
  *    0 \le z \perp Mz + q \ge 0
  * \f$
  * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ 
  * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. 
  * This sum is divided by \f$ \|q\| \f$ and then compared to tol.
  * \param[in] n size of the AVI
  * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
  * \param[in,out] error
  */
//  void avi_compute_error_only(int n,  double *z , double *w, double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

