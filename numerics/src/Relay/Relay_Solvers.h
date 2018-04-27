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
#ifndef RELAY_SOLVERS_H
#define RELAY_SOLVERS_H

/*!\file Relay_Solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon and Olivier Huber
  Subroutines for the resolution of relay problems.
*/

#include "RelayProblem.h"
#include "LinearComplementarityProblem.h"
#include "SolverOptions.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** set the default solver parameters and perform memory allocation for Relay
      \param[in] problem the RelayProblem structure which handles the problem (M,q)
      \param options the pointer to options to set
      \param solverId the identifier of the solver
  */
  int relay_setDefaultSolverOptions(RelayProblem* problem, SolverOptions* options, int solverId);


  /** relay_pgs is a projected Gauss-Seidel solver for relay problems.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param options the pointer to options to set
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Null diagonal term\n
   \author V. Acary
  */
  void relay_pgs(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for PGS
      \param options SolverOptions * the pointer to options to set
  */
  int relay_pgs_setDefaultSolverOptions(SolverOptions* options);

  /** relay_lexicolemke is a Lemke solver for  relay problems.\n
     * \param[in] problem structure that represents the Relay (M, q...)
     * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
     * \param options struct used to define the solver(s) and its (their) parameters
     * \param[out] info an integer which returns the termination value:\n
     0 = convergence,\n
     1 = no convergence,\n
     \author V. Acary
    */
  void relay_lexicolemke(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for Lemke
      \param options the pointer to options to set
  */
  int relay_lexicolemke_setDefaultSolverOptions(SolverOptions* options);

  /** relay_enum is enum solver for  relay problems.\n
     * \param[in] problem structure that represents the Relay (M, q...)
     * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
     * \param options struct used to define the solver(s) and its (their) parameters
     * \param[out] info an integer which returns the termination value:\n
     0 = convergence,\n
     1 = no convergence,\n
     2 = Null diagonal term\n
     \author V. Acary
    */
  void relay_enum(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for ENUM
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param options SolverOptions * the pointer to options to set
  */
  int relay_enum_setDefaultSolverOptions(RelayProblem* problem, SolverOptions* options);

  /** relay_path is a resolution of the Relay with its inherent MCP formulation and using path.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param options struct used to define the solver(s) and its (their) parameters
   * \param[out] info an integer which returns the termination value:\n
   *  0 = convergence,\n
   *  1 = no convergence,\n
   *  2 = Nul diagonal term\n
   * \author V. acary
   */
  void relay_path(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for PATH
      \param options the pointer to options to set
  */
  int relay_path_setDefaultSolverOptions(SolverOptions* options);

  /** Solve a Relay problem using the AVI framework and the solver by Cao
   * and Ferris.
   * \param[in] problem structure that represents the Relay (M, q, ...)
   * \param[in,out] z vector which on call is the initial point and on exit is the solution of the problem.
   * \param[in,out] w vector for computations
   * \param options struct used to define the solver(s) and its (their) parameters
   * \param[out] info an integer which returns the termination value:\n
   *  0 = convergence,\n
   *  1 = no convergence,\n
   * \author Olivier Huber
   */
  void relay_avi_caoferris(RelayProblem* problem, double* restrict z, double* restrict w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for
   * AVI_CAOFERRIS
   * \param options struct used to define the solver(s) and its (their) parameters
   */
  int relay_avi_caoferris_setDefaultSolverOptions(SolverOptions* options);

  /** Solve a Relay problem using the AVI framework and the solver by Cao and Ferris.
   * \warning This is only a test version. It does not take into account the
   * specificities of the problem like relay_avi_caoferris() does. Please do
   * not use this solver unless you have a pretty good reason.
   * \param[in] problem structure that represents the Relay (M, q, ...)
   * \param[in,out] z vector which on call is the initial point and on exit is the solution of the problem.
   * \param[in,out] w vector for computations
   * \param options struct used to define the solver(s) and its (their) parameters
   * \param[out] info an integer which returns the termination value:\n
   *  0 = convergence,\n
   *  1 = no convergence,\n
   * \author Olivier Huber
   */
  void relay_avi_caoferris_test(RelayProblem* problem, double* restrict z, double* restrict w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for AVI_CAOFERRIS
   * \param options struct used to define the solver(s) and its (their) parameters
   */
  int relay_avi_caoferris_test_setDefaultSolverOptions(SolverOptions* options);

  /** dr_latin is a specific latin (LArge Time INcrement)solver for dual relay problems.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Cholesky factorization failed,\n
   3 = Nul diagonal term\n
   * \param options struct used to define the solver(s) and its (their) parameters
   \author Nineb Sheherazade.
  */
  void dr_latin(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options)  ;

  /**  dr_nlgs is a specific nlgs (Non Linear Gauss Seidel) solver for dual relay problems.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Nul diagonal term\n
   * \param options struct used to define the solver(s) and its (their) parameters
  */
  void dr_nlgs(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution \n
     * of the LCP : \n
     * \f$
     *   -(Mz + q) \in  N_{[lb,ub]}(z)
     * \f$
     * The criterion is based on \f$ error = \|z- proj_{[lb,ub]}(z - \rho * (M*z+q)) \|, \rho >0\f$ \n
     * This error is divided by \f$ \|q\| \f$ and then compared to tol.
     * \param[in] problem structure that represents the Relay (M, q...)
     * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
     * \param[in] tolerance threshold used to validate the solution: if the error is less than this value, the solution is accepted
     * \param[in,out] error the actual error of the solution with respect to the problem
     * \return status: 0 : convergence, 1: error > tolerance
     * \author Vincent Acary
     */
  int relay_compute_error(RelayProblem* problem, double* restrict z , double* restrict w, double tolerance, double* restrict error);


  /** This function computes the projection on the boxr \f$ [lb,ub]\f$ of the vector  \f$z\f$  \n
     * \param[in,out] z a n-vector of doubles which returns the projection
     * \param[in,out] ub a n-vector of doubles which contains the upper bounds
     * \param[in,out] lb a n-vector of doubles which contains the lower bounds
     * \param[in,out] n size of the a n-vector
     * \author Vincent Acary
     */
  void project_on_box(int n, double* restrict z , double* restrict lb, double* restrict ub);

  /** This function transform a RelayProblem into a LinearComplementarityProblem
     * \param[in] problem A pointer to a Relay_problem to transform
     * \param[out] lcp_problem A pointer to a LinearComplementarity_problem resulting from the reformulation
     * \author Vincent Acary
     */
  void relay_to_lcp(RelayProblem* problem, LinearComplementarityProblem* lcp_problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
