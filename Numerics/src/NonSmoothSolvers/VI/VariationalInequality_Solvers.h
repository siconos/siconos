/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#ifndef VISOLVERS_H
#define VISOLVERS_H

/*!\file VariationalInequality_Solvers.h
  \brief Subroutines for the resolution of Variational Inequalites (VI) problems
*/

/*! \page VISolvers VI problems Solvers

This page gives an overview of the available solvers for Variational Inequality problems and their required parameters.

For each solver, the input argument are:
- a VariationalInequality
- the unknowns (x,fx)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions structure, which handles iparam and dparam
*/

#include "VariationalInequality.h"
#include "NumericsOptions.h"
#include "SolverOptions.h"
#include "VI_cst.h"
#include "SiconosCompat.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** General interface to solvers for variational inequality problem
  \param[in] problem the structure which handles the variational inequality problem
  \param[in,out] x global vector (n)
  \param[in,out] w global vector (n)
  \param[in,out] options structure used to define the solver(s) and their parameters
  \param[in] global_options for Numerics (verbose mode ...)
  \return result (0 if successful otherwise 1).
  */
  int variationalInequality_driver(VariationalInequality* problem, double *x , double *w, SolverOptions* options, NumericsOptions* global_options);

  /** set the default solver parameters and perform memory allocation for VariationalInequality
      \param options the pointer to the options to set
      \param solverId the identifier of the solver
  */
  int variationalInequality_setDefaultSolverOptions(SolverOptions* options, int solverId);

  /**Extra Gradient solver forvariational inequality problem based on the De Saxce Formulation
      \param problem the variational inequality problem to solve
      \param x global vector (n), in-out parameter
      \param w global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      dparam[3] : rho >0
  */
  void variationalInequality_ExtraGradient(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */
  
  int variationalInequality_ExtraGradient_setDefaultSolverOptions(SolverOptions* options);



  /** Fixed Point Projection solver for variational inequality problem based on the De Saxce Formulation
      \param problem the variational inequality problem to solve
      \param x global vector (n), in-out parameter
      \param w global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      dparam[3] : rho >0
  */
  void variationalInequality_FixedPointProjection(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options);


  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */
  int variationalInequality_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options);


  /** Hyperplane Projection solver for variational inequality problem based on the De Saxce Formulation
      \param problem the variational inequality problem to solve
      \param x global vector (n), in-out parameter
      \param w global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      dparam[3] : rho >0
  */
  void variationalInequality_HyperplaneProjection(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options);


  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */
  int variationalInequality_HyperplaneProjection_setDefaultSolverOptions(SolverOptions* options);


  /** VI Solver based on a merit function minimization with a line-search type
   * algorithm 
   * \param problem the variational inequality problem to solve
   * \param[in,out] x, as input, the initial guess; as output the solution if
   * the algorithm is successful
   * \param[in,out] F value of the function
   * \param info 0 if a solution is found
   * \param options the solver options
   */
  void variationalInequality_box_newton_QiLSA(VariationalInequality* problem, double *x, double *F, int* info, SolverOptions* options);

  /**  set the default solver parameters and perform memory allocation for a VI
   * solver
   * \param options the SolverOptions to set
   * \param solverId the id of the solver
   */
  int variationalInequality_common_setDefaultSolverOptions(SolverOptions* options, int solverId);



  /** Check for trivial solution in the variational inequality problem
      \param problem VariationalInequality*  the problem
      \param x global vector (n), in-out parameter
      \param fx global vector (n), in-out parameters
      \param options the pointer to the array of options to set
      \return info  =0 if a trivial solution has been found, else = -1
  */
  int checkTrivialCase_vi(VariationalInequality* problem , double* x, double* fx, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
