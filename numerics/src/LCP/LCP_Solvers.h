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
#ifndef LCP_SOLVERS_H
#define LCP_SOLVERS_H

/*!\file LCP_Solvers.h
  \brief Subroutines for the resolution of Linear Complementarity Problems.\n

  \author siconos-team@lists.gforge.inria.fr
*/

#include "NumericsFwd.h"
#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param[in] problem the LinearComplementarityProblem structure which handles the problem (M,q)
      \param options the pointer to the array of options to set
      \return info termination value
  */
  int linearComplementarity_setDefaultSolverOptions(LinearComplementarityProblem* problem, SolverOptions* options, int);

  /** lcp_qp uses a quadratic programm formulation for solving a LCP
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   *                0 : convergence  / minimization sucessfull\n
   *                1 : Too Many iterations\n
   *              2 : Accuracy insuficient to satisfy convergence criterion\n
   *                5 : Length of working array insufficient\n
   *                Other : The constraints are inconstent\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   * \author Vincent Acary
   */
  void lcp_qp(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_qp_setDefaultSolverOptions(SolverOptions* options);

  /** lcp_cpg is a CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0: convergence\n
   1: iter = itermax\n
   2: negative diagonal term\n
   3: pWp nul
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   \author Mathieu Renouf.
  */
  void lcp_cpg(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
    \param options the pointer to the array of options to set
  */
  int linearComplementarity_cpg_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for LCP.\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   \param[in,out] options structure used to define the solver and its parameters.
  */
  void lcp_pgs(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_pgs_setDefaultSolverOptions(SolverOptions* options);

  /** lcp_rpgs (Regularized Projected Gauss-Seidel ) is a solver for LCP, able to handle matrices with null diagonal terms.\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   * \param[in,out] options structure used to define the solver and its parameters.
   *

   \author Mathieu Renouf & Pascal Denoyelle
   \todo Sizing the regularization paramter and apply it only on null diagnal term

  */
  void lcp_rpgs(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_rpgs_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_psor Projected Succesive over relaxation solver for LCP. See cottle, Pang Stone Chap 5
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   * \param[in,out] options structure used to define the solver and its parameters.
   \todo use the relax parameter
   \todo add test
   \todo add a vector of relaxation parameter wtith an auto sizing (see SOR algorithm for linear solver.)

  */
  void lcp_psor(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_psor_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_nsqp use a quadratic programm formulation for solving an non symmetric LCP
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   *                0 : convergence  / minimization sucessfull\n
   *                1 : Too Many iterations\n
   *            2 : Accuracy insuficient to satisfy convergence criterion\n
   *                5 : Length of working array insufficient\n
   *                Other : The constraints are inconstent\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   * \author Vincent Acary
   *
   */
  void lcp_nsqp(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_nsqp_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_latin (LArge Time INcrements) is a basic latin solver for LCP.
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : Cholesky Factorization failed \n
   3 : nul diagonal term\n
   * \param[in,out] options structure used to define the solver and its parameters.

   \author Nineb Sheherazade.
  */
  void lcp_latin(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);
  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options  the pointer to the array of options to set
  */
  int linearComplementarity_latin_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_latin_w (LArge Time INcrements) is a basic latin solver with relaxation for LCP.\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : Cholesky Factorization failed \n
   3 : nul diagonal term\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *

   \author Nineb Sheherazade.
  */
  void lcp_latin_w(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_latin_w_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_lexicolemke is a direct solver for LCP based on pivoting method principle for degenerate problem \n
   * Choice of pivot variable is performed via lexicographic ordering \n
   *  Ref: "The Linear Complementarity Problem" Cottle, Pang, Stone (1992)\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   * 0 : convergence\n
   * 1 : iter = itermax\n
   * 2 : negative diagonal term\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   *\author Mathieu Renouf
   */
  void lcp_lexicolemke(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_lexicolemke_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_newton_min uses a nonsmooth Newton method based on the min formulation  (or max formulation) of the LCP
      \f$
      0 \le z \perp w \ge 0 \Longrightarrow \min(w,\rho z)=0 \Longrightarrow w = \max(0,w - \rho z)
      \f$

      \f$
      H(z) = H(\left[ \begin{array}{c} z \\ w \end{array}\right])= \left[ \begin{array}{c} w-Mz-q \\ min(w,\rho z) \end{array}\right] =0\\
      \f$


      References: Alart & Curnier 1990, Pang 1990
      * \param[in] problem structure that represents the LCP (M, q...)
      * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
      * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
      * \param[out] info an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : Problem in resolution in DGESV\n
      *                0 : convergence  / minimization sucessfull\n
      *                1 : Too Many iterations\n
      *              2 : Accuracy insuficient to satisfy convergence criterion\n
      *                5 : Length of working array insufficient\n
      *                Other : The constraints are inconstent\n
      * \param[in,out] options structure used to define the solver and its parameters.
      *
      \author Vincent Acary

      \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
      \todo Add rules for the computation of the penalization rho
      \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
  */
  void lcp_newton_min(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_newton_min_setDefaultSolverOptions(SolverOptions* options);


  /** lcp_newton_FB use a nonsmooth newton method based on the Fischer-Bursmeister convex function
   *
   *
   * \f$
   *   0 \le z \perp w \ge 0 \Longrightarrow \phi(z,w)=\sqrt{z^2+w^2}-(z+w)=0
   * \f$

   * \f$
   *   \Phi(z) = \left[ \begin{array}{c}  \phi(z_1,w_1) \\ \phi(z_1,w_1) \\ \vdots \\  \phi(z_n,w_n)  \end{array}\right] =0\\
   * \f$
   *
   *
   * References: Alart & Curnier 1990, Pang 1990
   *
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *                2 - failure in the descent direction search (in LAPACK) \n
   *
   * \param[in,out] options structure used to define the solver and its parameters.
   * \author Vincent Acary and Olivier Huber
   *
   * \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
   * \todo Add rules for the computation of the penalization rho
   * \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
   */

  void lcp_newton_FB(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** lcp_newton_minFB use a nonsmooth newton method based on both a min and Fischer-Bursmeister reformulation
   * References: Facchinei--Pang (2003)
   *
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *                2 - failure in the descent direction search (in LAPACK) \n
   *
   * \param[in,out] options structure used to define the solver and its parameters.
   * \author Olivier Huber
   */
  void lcp_newton_minFB(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /**
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   0 : convergence\n
   1 : iter = itermax\n
   2 : negative diagonal term
   * \param[in,out] options structure used to define the solver and its parameters.

   \author Olivier Bonnefon
  */
  void lcp_path(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** enumerative solver
  * \param[in] problem structure that represents the LCP (M, q...)
  * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
  * \param[out] info an integer which returns the termination value:\n
  0 : success\n
  1 : failed\n
  * \param[in,out] options structure used to define the solver and its parameters.

  \author Olivier Bonnefon
  */
  void lcp_enum(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /**
  * \param[in] problem structure that represents the LCP (M, q...)
  * \param[in,out] options structure used to define the solver and its parameters.
  * \param[in] withMemAlloc If it is not 0, then the necessary work memory is allocated.
  */
  void lcp_enum_init(LinearComplementarityProblem* problem, SolverOptions* options, int withMemAlloc);
  /**
  * \param[in] problem structure that represents the LCP (M, q...)
  * \param[in,out] options structure used to define the solver and its parameters.
  * \param[in]  withMemAlloc If it  is not 0, then the work memory is free.
  */
  void lcp_enum_reset(LinearComplementarityProblem* problem, SolverOptions* options, int withMemAlloc);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param problem structure that represents the LCP (M, q...)
      \param  options  the pointer to the array of options to set
  */
  int linearComplementarity_enum_setDefaultSolverOptions(LinearComplementarityProblem* problem, SolverOptions* options);

  /** lcp_avi_caoferris is a direct solver for LCP based on an Affine Variational Inequalities (AVI) reformulation\n
   * The AVI solver is here the one from Cao and Ferris \n
   *  Ref: "A Pivotal Method for Affine Variational Inequalities" Menglin Cao et Michael Ferris (1996)\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   * 0 : convergence\n
   * 1 : iter = itermax\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   *\author Olivier Huber
   */
  void lcp_avi_caoferris(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** lcp_pivot is a direct solver for LCP based on a pivoting method\n
   * It can currently use Bard, Murty's least-index or Lemke rule for choosing
   * the pivot. The default one is Lemke and it cam be changed by setting
   * iparam[2]. The list of choices are in the enum LCP_PIVOT (see lcp_cst.h).
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   * 0 : convergence\n
   * 1 : iter = itermax\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   *\author Olivier Huber
   */
  void lcp_pivot(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);
  void lcp_pivot_covering_vector(LinearComplementarityProblem* problem, double* u , double* s, int *info , SolverOptions* options, double* cov_vec);
  void lcp_pivot_lumod(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);
  void lcp_pivot_lumod_covering_vector(LinearComplementarityProblem* problem, double* u , double* s, int *info , SolverOptions* options, double* cov_vec);

  /** lcp_pathsearch is a direct solver for LCP based on the pathsearch algorithm\n
   * \warning this solver is available for testing purposes only! consider
   * using lcp_pivot() if you are looking for simular solvers
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   * 0 : convergence\n
   * 1 : iter = itermax\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   *\author Olivier Huber
   */
  void lcp_pathsearch(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** lcp_gams uses the solver provided by GAMS \n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   * 0 : convergence\n
   * 1 : iter = itermax\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   *\author Olivier Huber
   */
  void lcp_gams(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** generic interface used to call any LCP solver applied on a Sparse-Block structured matrix M, with a Gauss-Seidel process
   * to solve the global problem (formulation/solving of local problems for each row of blocks)
   * \param[in] problem structure that represents the LCP (M, q...). M must be a SparseBlockStructuredMatrix
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param info an integer which returns the termination value:\n
   0 : convergence\n
   >0 : failed, depends on local solver
   * \param[in,out] options structure used to define the solver and its parameters.
   * \author Mathieu Renouf, Pascal Denoyelle, Franck Perignon
   */
  void lcp_nsgs_SBM(LinearComplementarityProblem* problem, double *z, double *w, int* info, SolverOptions* options);
  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to  the array of options to set
  */
  int linearComplementarity_nsgs_SBM_setDefaultSolverOptions(SolverOptions* options);

  /** Construct local problem from a "global" one
   * \param rowNumber index of the local problem
   * \param blmat matrix containing the problem
   * \param local_problem problem to fill
   * \param q big q
   * \param z big z
   */
  void lcp_nsgs_SBM_buildLocalProblem(int rowNumber, const SparseBlockStructuredMatrix* const blmat, LinearComplementarityProblem* local_problem, double* q, double* z);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution \n
   * of the LCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[in] tolerance threshold used to validate the solution: if the error is less than this value, the solution is accepted
   * \param[out] error the actual error of the solution with respect to the problem
   * \return status: 0 : convergence, 1: error > tolerance
   * \author Pascal Denoyelle, Franck Perignon
   */
  int lcp_compute_error(LinearComplementarityProblem* problem, double *z , double *w, double tolerance, double* error);

  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution \n
  * of the LCP : \n
  * \f$
  *    0 \le z \perp Mz + q \ge 0
  * \f$
  * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
  * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
  * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
  * \param[in] n size of the LCP
  * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
  * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
  * \param[out] error the result of the computation
  * \author Pascal Denoyelle, Franck Perignon
  */
  void lcp_compute_error_only(unsigned int n,  double *z , double *w, double * error);

  /*   /\** Function used to extract from LCP matrix the part which corresponds to non null z */
  /*    *\/ */
  /*   int extractLCP( NumericsMatrix* MGlobal, double *z , int *indic, int *indicop, double *submatlcp , double *submatlcpop, */
  /*     int *ipiv , int *sizesublcp , int *sizesublcpop); */

  /*   /\** Function used to solve the problem with sub-matrices from extractLCP */
  /*    *\/ */
  /*   int predictLCP(int sizeLCP, double* q, double *z , double *w , double tol, */
  /*   int *indic , int *indicop , double *submatlcp , double *submatlcpop , */
  /*    int *ipiv , int *sizesublcp , int *sizesublcpop , double *subq , double *bufz , double *newz); */

  /** Interface to solvers for Linear Complementarity Problems, dedicated to dense matrix storage
      \param[in] problem the LinearComplementarityProblem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
      \author Nineb Sheherazade, Mathieu Renouf, Franck Perignon
  */
  int lcp_driver_DenseMatrix(LinearComplementarityProblem* problem, double *z , double *w, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

