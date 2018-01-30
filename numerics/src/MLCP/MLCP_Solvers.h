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
- a MixedLinearComplementarityProblem
- the unknowns (z,w)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions structure, which handles iparam and dparam

\section mlcpPGS PGS Solver
Projected Gauss-Seidel solver

 function: mlcp_pgs() \n
 parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- iparam[2] (in): 0 for implicit, 1 for explicit
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error

\section mlcpPGS_SBM PGS Solver for SBM storage
Projected Gauss-Seidel solver
- iparam[0] (in): maximum number of iterations allowed for GS process
- iparam[1] (out): number of GS iterations processed
- iparam[2] (out): sum of all local number of iterations (if it has sense for the local solver)
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): sum of all local error values

\section mlcpRPGS RPGS Solver
Regularized Projected Gauss-Seidel, solver for MLCP, able to handle with matrices with null diagonal terms

 function: mlcp_rpgs() \n
 parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): rho

\section mlcpPSOR PSOR Solver
Projected Succesive over relaxation solver for MLCP. See cottle, Pang Stone Chap 5
 function: mlcp_psor() \n
 parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): omega

\section mlcpRPSOR RPSOR Solver
Regularized Projected Succesive over relaxation solver for MLCP
 function: mlcp_rpsor() \n
 parameters:
- iparam[0] (in): maximum number of iterations allowed
- iparam[1] (out): number of iterations processed
- dparam[0] (in): tolerance
- dparam[1] (out): resulting error
- dparam[2] (in): omega
- dparam[3] (in): rho

\section mlcpPath Path (Ferris) Solver
The path solver must be initialize:\n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n
 function: mlcp_path() \n
 parameters:
- dparam[0] (in): tolerance

\section mlcpENUM ENUM Solver
The enumeratif solver must be initialize:\n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n

 function: mlcp_enum() \n
The enumeratif solver must be initialize: \n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n
 parameters:
- dparam[0] (in): a positive value, tolerane about the sign.
- iparam[4] (in) :  use DGELS (1) or DGESV (0).
- dWork : working float zone size : The number of doubles is retruned by the function  mlcp_driver_get_dwork() . MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function  mlcp_driver_get_iwork() . MUST BE ALLOCATED BY THE USER.

\section mlcpSIMPLEX SIMPLEX Solver
The simplex solver must be initialize:\n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n

 function: mlcp_simplex() \n
The simplex solver must be initialize: \n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n

 parameters:
- iparam[0] (in): Max number of iteration (example: 1000000).
- dparam[0] (in): A positive value, tolerance to consider that a var is null(ex: 10e-12).
- dparam[1] (in): A positive value, tolerance to consider that complementarity holds(ex: 10e-12).
- dparam[2] (in): A positive value, tolerance to consider that a var is negative(ex: 10e-9).
- dWork : working float zone size : The number of doubles is retruned by the function  mlcp_driver_get_dwork() . MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function  mlcp_driver_get_iwork() . MUST BE ALLOCATED BY THE USER.

\section mlcpDIRECT_ENUM DIRECT_ENUM Solver
The direct and enumeratif solver must be initialize: \n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n

 function: mlcp_direct_enum() \n
 parameters:
- iparam[5] (in): Number of registered configurations.
- iparam[7] (out): Number of case the direct solved failed.
- dparam[0] (in): A positive value, tolerane about the sign.
- dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).
- dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).
- dWork : working float zone size : The number of doubles is retruned by the function  mlcp_driver_get_dwork() . MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function  mlcp_driver_get_iwork() . MUST BE ALLOCATED BY THE USER.

\section mlcpDIRECT_PATH DIRECT_PATH Solver
The path solver must be initialize: \n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n

 function: mlcp_direct_path() \n
 parameters:
- iparam[5] (in): Number of registered configurations.
- iparam[7] (out): Number of case the direct solved failed.
- dparam[0] (in): Tolerance.
- dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).
- dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).
- dWork : working float zone size : The number of doubles is retruned by the function  mlcp_driver_get_dwork() . MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function  mlcp_driver_get_iwork() . MUST BE ALLOCATED BY THE USER.

\section mlcpDIRECT_SIMPLEX DIRECT_SIMPLEX Solver
The direct and simplex solver must be initialize: \n
1) Initialize the solver with  mlcp_driver_init() .\n
2) Use a lot with  mlcp_driver() .\n
3) Reset the solver with  mlcp_driver_reset() .\n
 function: mlcp_direct_simplex() \n
 parameters:
- iparam[0] (in): Max number of iteration (example: 1000000).
- iparam[5] (in): Number of registered configurations.
- iparam[7] (out): Number of case the direct solved failed.
- dparam[0] (in): A positive value, tolerance to consider that a var is null(ex: 10e-12).
- dparam[1] (in): A positive value, tolerance to consider that complementarity holds(ex: 10e-12).
- dparam[2] (in): A positive value, tolerance to consider that a var is negative(ex: 10e-9).
- dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).
- dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).
- dWork : working float zone size : The number of doubles is retruned by the function  mlcp_driver_get_dwork() . MUST BE ALLOCATED BY THE USER.
- iWork : working int zone size : . The number of double is retruned by the function  mlcp_driver_get_iwork() . MUST BE ALLOCATED BY THE USER.

*/

#include "SolverOptions.h"
#include "MixedLinearComplementarityProblem.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

#include "mlcp_enum.h"
#include "mlcp_simplex.h"
#include "mlcp_path_enum.h"
#include "mlcp_direct_path_enum.h"
#include "mlcp_direct_enum.h"
#include "mlcp_direct_simplex.h"
#include "mlcp_direct_path.h"
#include "mlcp_direct_FB.h"
#include "mlcp_FB.h"

#include "mlcp_cst.h"


  /** General interface to initialize a solver.\n
      Must be call for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \author Olivier Bonnefon
  */
  void mlcp_driver_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
  /** General interface to reset a solver.\n
      Must be call for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \author Olivier Bonnefon
  */
  void mlcp_driver_reset(MixedLinearComplementarityProblem* problem, SolverOptions* options);
  /*
    Alloc the needed working memory if it was not yet be done, ie options->iwork and options->dwork must be null, else do nothing.
    \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
    \param[in] options structure used to define the solver(s) and their parameters
    \return 0 iff the memory has not be done.
   */
  int mlcp_alloc_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options);
  /*
    Free the working memory.
    \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
    \param[in] options structure used to define the solver(s) and their parameters
   */
  void mlcp_free_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options);
  /** General interface to get the number of integers that must be allocated for the solver.\n
      Must be use for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \return the number of integers that must be allocated by the user.

      \author Olivier Bonnefon
  */
  int mlcp_driver_get_iwork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
  /** General interface to get the number of doubles that must be allocated for the solver.\n
      Must be use for the following solvers:\n
      - mlcp_enum
      - mlcp_path
      - mlcp_simplex
      - mlcp_direct_enum
      - mlcp_direct_path
      - mlcp_direct_simplex

      \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
      \param[in] options structure used to define the solver(s) and their parameters
      \return the number of doubles that must be allocated by the user.

      \author Olivier Bonnefon
  */
  int mlcp_driver_get_dwork(MixedLinearComplementarityProblem* problem, SolverOptions* options);

  /**  mlcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in,out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in,out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_pgs(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /**  mlcp_rpgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in,out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in,out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_rpgs(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** mlcp_psor (projected successive overrelaxation method) is a solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in,out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in,out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_psor(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** mlcp_rpsor (regularized projected successive overrelaxation method) is a solver for MLCP.
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in,out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in,out] options structure used to define the solver and its parameters.
   \author Vincent Acary
  */
  void mlcp_rpsor(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** path solver
   * \param[in] problem structure that represents the MLCP (n,m,M, q...)
   * \param[in,out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in,out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_path(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** enum solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in,out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** direct solver
   * \param[in]  problem MixedLinearComplementarityProblem* problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in,out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_direct(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** direct-enum solver
  * \param[in] problem structure that represents the MLCP (n,mM, q...)
  * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
  * \param[out] info an integer which returns the termination value:\n
  0 : success,it found a solution\n
  1 : echec,it did not find any solution\n
  \param[in,out] options structure used to define the solver and its parameters.
  \author Olivier Bonnefon
  */
  void mlcp_direct_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** direct-simplex solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in,out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_direct_simplex(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** direct-path solver
  * \param[in] problem structure that represents the MLCP (n,mM, q...)
  * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
  * \param[out] info an integer which returns the termination value:\n
  0 : success,it found a solution\n
  1 : echec,it did not find any solution\n
  \param[in,out] options structure used to define the solver and its parameters.
  \author Olivier Bonnefon
  */
  void mlcp_direct_path(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);


  /** simplex solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in,out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_simplex(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** Fischer Burmeister solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in,out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_FB(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);
  /** Direct Fischer Burmeister solver
   * \param[in] problem structure that represents the MLCP (n,mM, q...)
   * \param[out] z a m+n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[out] w a m+n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : success,it found a solution\n
   1 : echec,it did not find any solution\n
   \param[in,out] options structure used to define the solver and its parameters.
   \author Olivier Bonnefon
  */
  void mlcp_direct_FB(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);


  /** generic interface used to call any MLCP solver applied on a Sparse-Block structured matrix M, with a Gauss-Seidel process
   * to solve the global problem (formulation/solving of local problems for each row of blocks)
   * \param[in] problem structure that represents the LCP (M, q...). M must be a SparseBlockStructuredMatrix
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param info an integer which returns the termination value:\n
   0 : convergence\n
   >0 : failed, depends on local solver
   * \param[in,out] options structure used to define the solver and its parameters.
   * \author Vincent Acary
   */
  void mlcp_pgs_SBM(MixedLinearComplementarityProblem* problem, double *z, double *w, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param problem structure that represents the MLCP (M, q...). M must be a SparseBlockStructuredMatrix
      \param options the pointer to  the array of options to set
  */
  int mixedLinearComplementarity_pgs_SBM_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* options);

  // need a svn add mlcp_GaussSeidel_SBM ...
  //  void mlcp_GaussSeidel_SBM(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options, int numberOfSolvers);

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
     \param[in] tolerance threshold used to validate the solution: if the error is less than this value, the solution is accepted
     \param[in,out] error
     \return status: 0 : convergence, 1: error > tolerance
    \author Vincent Acary form the routine  filter_result_LCP.c of Pascal Denoyelle
  */
  int mlcp_compute_error(MixedLinearComplementarityProblem* problem, double *z, double *w, double tolerance, double * error);

  /*
    Default memory allocator.
      \param Problem * the pointer to the array of options to set.
      \param SolverOptions * the pointer to option.
  */
  void  mixedLinearComplementarity_default_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions);
  /*
    Default memory manager to free the memory located in the options.
      \param problem * the pointer to the array of options to set.
      \param SolverOptions * the pointer to option.
   */
  void  mixedLinearComplementarity_deleteDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions);

  /** set the default solver parameters and perform memory allocation for MixedLinearComplementarity
      \param problem * the pointer to the array of options to set.
      \param pOptions SolverOptions * the pointer to option.
  */
  int mixedLinearComplementarity_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions);
  int mixedLinearComplementarity_directEnum_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_directFB_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_directPath_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_directPathEnum_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_directSimplex_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_enum_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_fb_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_path_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_pathEnum_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_pgs_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_rpgs_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_simplex_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_rpsor_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);
  int mixedLinearComplementarity_psor_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
