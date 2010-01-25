/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#ifndef PRSOLVERS_H
#define PRSOLVERS_H

/*!\file Relay_Solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
  Subroutines for the resolution of relay problems.
*/

/*! \page RelaySolvers Relay Problems Solvers

This page gives an overview of the available solvers for relay problems and their required parameters.

For each solver, the input argument are:
- a Relay_Problem
- the unknowns (z,w)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a Solver_Options structure, which handles iparam and dparam

Note: function names starting with dr -> dual relay, with pr: primal relay

\section relayLatin Latin
LArge Time INcrements solver

\bf function: pr_latin() or dr_latin() \n
\bf parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): latin parameter

\section relayNLGS Non-linear Gauss Seidel
LArge Time INcrements solver

\bf function: pr_nlgs() or dr_nlgs()\n
\bf parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error

*/

#include "Relay_Problem.h"
#include "LinearComplementarity_Problem.h"
#include "Solver_Options.h"

#ifdef __cplusplus
extern "C" {
#endif


  /** General interface to solver for primal-relay problems
      \param[in] problem the Relay_Problem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
   * \author Nineb Sheherazade.
   */
  int relay_driver(Relay_Problem* problem, double *z , double *w, Solver_Options* options,  Numerics_Options* global_options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param Solver_Options * the pointer to options to set
      \param char * the string which identify the solver
  */
  int relay_setDefaultSolverOptions(Relay_Problem* problem, Solver_Options* options, char *);


  /** relay_nlgs is a projected Gauss-Seidel solver for relay problems.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Null diagonal term\n
   \author V. Acary
  */
  void relay_pgs(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** set the default solver parameters and perform memory allocation for PGS
      \param Solver_Options * the pointer to options to set
  */
  int relay_pgs_setDefaultSolverOptions(Solver_Options* options);

  /** relay_lexicolemke is a Lemke solver for  relay problems.\n
     * \param[in] problem structure that represents the Relay (M, q...)
     * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
     * \param[out] info an integer which returns the termination value:\n
     0 = convergence,\n
     1 = no convergence,\n
     2 = Null diagonal term\n
     \author V. Acary
    */
  void relay_lexicolemke(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options, Numerics_Options* global_options);

  /** set the default solver parameters and perform memory allocation for Lemke
      \param Solver_Options * the pointer to options to set
  */
  int relay_lexicolemke_setDefaultSolverOptions(Solver_Options* options);

  /** relay_nlgs is enum solver for  relay problems.\n
     * \param[in] problem structure that represents the Relay (M, q...)
     * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
     * \param[out] info an integer which returns the termination value:\n
     0 = convergence,\n
     1 = no convergence,\n
     2 = Null diagonal term\n
     \author V. Acary
    */
  void relay_enum(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options, Numerics_Options* global_options);

  /** set the default solver parameters and perform memory allocation for ENUM
      \param Solver_Options * the pointer to options to set
  */
  int relay_enum_setDefaultSolverOptions(Relay_Problem* problem, Solver_Options* options);

  /** relay_path is a resolution of the Relay with its inherent MCP formulation and using path.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Nul diagonal term\n
   \author V. acary
  */
  void relay_path(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options);


  /** pr_latin is a specific latin solver for primal relay problems.
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Cholesky factorization failed,\n
   3 = Nul diagonal term\n
   \author Nineb Sheherazade.
  */
  void pr_latin(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options);


  /** dr_latin is a specific latin (LArge Time INcrement)solver for dual relay problems.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Cholesky factorization failed,\n
   3 = Nul diagonal term\n
   \author Nineb Sheherazade.
  */
  void dr_latin(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options)  ;

  /**  dr_nlgs is a specific nlgs (Non Linear Gauss Seidel) solver for dual relay problems.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Nul diagonal term\n
   \author Nineb Sheherazade.
  */
  void dr_nlgs(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** pr_gsnl is a specific gsnl (Gauss Seidel Non Linear)solver for relay problems.
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Nul diagonal term\n
   \author Nineb Sheherazade.
  */
  void pr_gsnl(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution \n
     * of the LCP : \n
     * \f$
     *   -(Mz + q) \in  N_{[lb,ub]}(z)
     * \f$
     * The criterion is based on \f$ error = \|z- proj_{[lb,ub]}(z - \rho * (M*z+q)) \|, \rho >0\f$ \n
     * This error is divided by \f$ \|q\| \f$ and then compared to tol.\n
     * \param[in] problem structure that represents the Relay (M, q...)
     * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
     * \param[in] tolerance
     * \param[in,out] error
     * \return status: 0 : convergence, 1: error > tolerance
     * \author Vincent Acary
     */
  int relay_compute_error(Relay_Problem* problem, double *z , double *w, double tolerance, double* error);


  /** This function computes the projection on the boxr \f$ [lb,ub]\f$ of the vector  \f$z\f$  \n
     * \param[in,out] z a n-vector of doubles which returns the projection
     * \param[in,out] ub a n-vector of doubles which contains the upper bounds
     * \param[in,out] lb a n-vector of doubles which contains the lower bounds
     * \return status: 0
     * \author Vincent Acary
     */
  int projectiononbox(double *z , double *lb, double * ub, int n);

  /** This function transform a Relay_Problem into a LinearComplementarity_Problem
     * \param[in] A pointer to a Relay_problem to transform
     * \param[out] A pointer to a LinearComplementarity_problem resulting from the reformulation
     * \author Vincent Acary
     */


  void relay_tolcp(Relay_Problem* problem, LinearComplementarity_Problem * lcp_problem);

#ifdef __cplusplus
}
#endif

#endif
