/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
  \author Vincent Acary \n
  Last Modifications : Olivier Bonnefon

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
- iparam[2] (in): 0 for implicit, 1 for explicit
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
The path solver must be initialize:\n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n
\bf function: mlcp_path() \n
\bf parameters:
- dparam[0] (in): tolerance

\section mlcpENUM ENUM Solver
The enumeratif solver must be initialize:\n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n

\bf function: mlcp_enum() \n
The enumeratif solver must be initialize: \n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n
\bf parameters:
- dparam[0] (in): a positive value, tolerane about the sign.
- dWork : working float zone size : The number of doubles is retruned by the function \bf mlcp_driver_get_dwork() \bf. MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function \bf mlcp_driver_get_iwork() \bf. MUST BE ALLOCATED BY THE USER.

\section mlcpSIMPLEX SIMPLEX Solver
The simplex solver must be initialize:\n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n

\bf function: mlcp_simplex() \n
The simplex solver must be initialize: \n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n

\bf parameters:
- iparam[0] (in): Max number of iteration (example: 1000000).
- dparam[0] (in): A positive value, tolerance to consider that a var is null(ex: 10e-12).
- dparam[1] (in): A positive value, tolerance to consider that complementarity holds(ex: 10e-12).
- dparam[2] (in): A positive value, tolerance to consider that a var is negative(ex: 10e-9).
- dWork : working float zone size : The number of doubles is retruned by the function \bf mlcp_driver_get_dwork() \bf. MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function \bf mlcp_driver_get_iwork() \bf. MUST BE ALLOCATED BY THE USER.

\section mlcpDIRECT_ENUM DIRECT_ENUM Solver
The direct and enumeratif solver must be initialize: \n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n

\bf function: mlcp_direct_enum() \n
\bf parameters:
- iparam[5] (in): Number of registered configurations.
- iparam[7] (out): Number of case the direct solved failed.
- dparam[0] (in): A positive value, tolerane about the sign.
- dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).
- dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).
- dWork : working float zone size : The number of doubles is retruned by the function \bf mlcp_driver_get_dwork() \bf. MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function \bf mlcp_driver_get_iwork() \bf. MUST BE ALLOCATED BY THE USER.

\section mlcpDIRECT_PATH DIRECT_PATH Solver
The path solver must be initialize: \n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n

\bf function: mlcp_direct_path() \n
\bf parameters:
- iparam[5] (in): Number of registered configurations.
- iparam[7] (out): Number of case the direct solved failed.
- dparam[0] (in): Tolerance.
- dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).
- dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).
- dWork : working float zone size : The number of doubles is retruned by the function \bf mlcp_driver_get_dwork() \bf. MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function \bf mlcp_driver_get_iwork() \bf. MUST BE ALLOCATED BY THE USER.

\section mlcpDIRECT_SIMPLEX DIRECT_SIMPLEX Solver
The direct and simplex solver must be initialize: \n
1) Initialize the solver with \bf mlcp_driver_init() \bf.\n
2) Use a lot with \bf mlcp_driver() \bf.\n
3) Reset the solver with \bf mlcp_driver_reset() \bf.\n
\bf function: mlcp_direct_simplex() \n
\bf parameters:
- iparam[0] (in): Max number of iteration (example: 1000000).
- iparam[5] (in): Number of registered configurations.
- iparam[7] (out): Number of case the direct solved failed.
- dparam[0] (in): A positive value, tolerance to consider that a var is null(ex: 10e-12).
- dparam[1] (in): A positive value, tolerance to consider that complementarity holds(ex: 10e-12).
- dparam[2] (in): A positive value, tolerance to consider that a var is negative(ex: 10e-9).
- dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).
- dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).
- dWork : working float zone size : The number of doubles is retruned by the function \bf mlcp_driver_get_dwork() \bf. MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function \bf mlcp_driver_get_iwork() \bf. MUST BE ALLOCATED BY THE USER.

*/

#include "Numerics_Options.h"
#include "Solver_Options.h"
#include "MixedLinearComplementarity_Problem.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "mlcp_enum.h"
#include "mlcp_simplex.h"
#include "mlcp_direct_enum.h"
#include "mlcp_direct_simplex.h"
#include "mlcp_direct_path.h"
#include "mlcp_direct_FB.h"
#include "mlcp_FB.h"

  /** General interface to initialize a solver.\n
      Must be call for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarity_Problem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \author Olivier Bonnefon
  */
  void mlcp_driver_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
  /** General interface to reset a solver.\n
      Must be call for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarity_Problem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \author Olivier Bonnefon
  */
  void mlcp_driver_reset(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

  /** General interface to get the number of integers that must be allocated for the solver.\n
      Must be use for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarity_Problem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \return the number of integers that must be allocated by the user.

      \author Olivier Bonnefon
  */
  int mlcp_driver_get_iwork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
  /** General interface to get the number of doubles that must be allocated for the solver.\n
      Must be use for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarity_Problem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \return the number of doubles that must be allocated by the user.

      \author Olivier Bonnefon
  */
  int mlcp_driver_get_dwork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

  /**  mlcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in-out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_pgs(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /**  mlcp_rpgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in-out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_rpgs(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** mlcp_psor (projected successive overrelaxation method) is a solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in-out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_psor(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** mlcp_rpsor (regularized projected successive overrelaxation method) is a solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in-out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_rpsor(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** path solver
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in-out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in-out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_path(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** enum solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_enum(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** direct solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_direct(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** direct-enum solver
  * \param[in] problem structure that represents the MLCP (n,mM, q...)
  * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
  * \param[out] info an integer which returns the termination value:\n
  0 : success,it found a solution\n
  1 : echec,it did not find any solution\n
  \param[in-out] options structure used to define the solver and its parameters.
  \author Olivier Bonnefon
  */
  void mlcp_direct_enum(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** direct-simplex solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_direct_simplex(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** direct-path solver
  * \param[in] problem structure that represents the MLCP (n,mM, q...)
  * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
  * \param[out] info an integer which returns the termination value:\n
  0 : success,it found a solution\n
  1 : echec,it did not find any solution\n
  \param[in-out] options structure used to define the solver and its parameters.
  \author Olivier Bonnefon
  */
  void mlcp_direct_path(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);


  /** simplex solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_simplex(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /** Fischer Burmeister solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_FB(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);
  /** Direct Fischer Burmeister solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in-out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_direct_FB(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);

  /**
    This function checks the validity of the vector z as a solution \n
    of the MLCP : \n
    \f{eqnarray*}
     \left\lbrace
      \begin{array}{l}
      A u + Cv +a =w1\\
      D u + Bv +b = w2\\
      0 \le v \perp  w2 \ge 0\\
      w1=0\\
      \end{array}
     \right.
     w =
     \left\lbrace
     \begin{array}{c}
     w1\\
     w2\\
     \end{array}
     \right\rbrace
     \f}
     The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
     with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
     This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
     It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
     \param[in] problem structure that represents the MLCP (n,m,M, q... or (A,B,C...))
     \param[in,out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
     \param[in,out] w a m+n-vector of doubles which returns the solution of the problem.
     \param[in] tolerance
     \param[in,out] error
     \return status: 0 : convergence, 1: error > tolerance
    \author Vincent Acary form the routine  filter_result_LCP.c of Pascal Denoyelle
  */
  int mlcp_compute_error(MixedLinearComplementarity_Problem* problem, double *z, double *w, double tolerance, double * error);




#ifdef __cplusplus
}
#endif

#endif
