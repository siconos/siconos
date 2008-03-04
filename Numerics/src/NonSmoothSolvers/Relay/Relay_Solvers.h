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
#include "Solver_Options.h"

#ifdef __cplusplus
extern "C" {
#endif

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

  /** pr_nlgs is a specific nlgs(non linear Gauss-Seidel) solver for primal relay problems.\n
   * \param[in] problem structure that represents the Relay (M, q...)
   * \param[in-out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 = convergence,\n
   1 = no convergence,\n
   2 = Nul diagonal term\n
   \author Nineb Sheherazade.
  */
  void pr_nlgs(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

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

#ifdef __cplusplus
}
#endif

#endif
