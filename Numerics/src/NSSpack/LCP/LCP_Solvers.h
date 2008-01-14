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
#ifndef LCP_SOLVERS_H
#define LCP_SOLVERS_H

/*! \page LCPSolvers Linear Complementarity problems (LCP)
  \section lcpIntro The problem
  Find \f$(z,w)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ z + q = w \\
  0 \le w \perp z \ge 0 \\
  \end{array}
  \right.
  \f$

  \f$ w, z, q\f$ are vectors of size n and \f$ M \f$ is a nXn matrix.

  \section lcpSolversList Available solvers
  Use the generic function lcp_solver() or lcp_solver_block(), to call one the the specific solvers listed below:

  - lcp_qp(), quadratic programm formulation
  - lcp_cpg(), CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.
  - lcp_pgs(), PGS is a basic Projected Gauss-Seidel solver for LCP.
  - lcp_rpgs(), Regularized Projected Gauss-Seidel, is a solver for LCP, able to handle matrices with null diagonal terms
  - lcp_psor(), Projected Succesive over relaxation solver for LCP. See cottle, Pang Stone Chap 5
  - lcp_nsqp(), quadratic programm formulation for solving an non symmetric LCP
  - lcp_latin(), (LArge Time INcrements) is a basic latin solver for LCP.
  - lcp_latin_w(), (LArge Time INcrements) is a basic latin solver with relaxation for LCP
  - lcp_lexicolemke(), direct solver for LCP based on pivoting method principle for degenerate problem.\n
  Choice of pivot variable is performed via lexicographic ordering
  Ref: "The Linear Complementary Problem" Cottle, Pang, Stone (1992)\n
  - lcp_newton_min(), nonsmooth Newton method based on the min formulation  (or max formulation) of the LCP
  - lcp_newton_FB(), uses a nonsmooth newton method based on the Fischer-Bursmeister convex function

  The structure method, argument of lcp_solver(), is used to give the name and parameters of the required solver.

  (see the functions/solvers list in lcp_solvers.h)

*/

/*!\file LCP_Solvers.h
  Subroutines for the resolution of LCP problems.\n

  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

#include "NSSTools.h"

/** specific method for LCP
    \param name       name of the solver.
    \param itermax    maximum number of iterations.
    \param tol        convergence criteria value.
    \param k_latin    latin coefficient.
    \param relax      relaxation coefficient.
    \param rho        regularization coefficient
    \param chat       output boolean ( 0 = no output log ).
    \param normType   name norm (not yet available).
    \param iter       final number of iterations.
    \param err        final value of error criteria.
*/
typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  double relax;
  double rho;
  int    chat;
  char   normType[64];
  int    iter;
  double err;

} method_lcp;

#ifdef __cplusplus
extern "C" {
#endif


  /** lcp_qp uses a quadratic programm formulation for solving a LCP
   * \param nn    On enter, an integer which represents the dimension of the system.
   * \param M     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
   * \param q     On enter, a nn-vector of doubles which contains the components of the right hand side vector.
   * \param z     On return, a nn-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param w     On return, a nn-vector of doubles which returns the solution of the problem.
   * \param info  On return, an integer which returns the termination value:\n
   *                0 : convergence  / minimization sucessfull\n
   *                1 : Too Many iterations\n
   *              2 : Accuracy insuficient to satisfy convergence criterion\n
   *                5 : Length of working array insufficient\n
   *                Other : The constraints are inconstent\n
   *
   * \param iparamLCP    On enter/return, vector of integers (not used here).
   * \param dparamLCP    On enter/return, vector of doubles:\n
   *                - dparamLCP[0] = tol    On enter, the tolerance required.
   *
   * \author Vincent Acary
   */
  void lcp_qp(int *nn , double *M , double *q , double *z , double *w , int *info ,
              int *iparamLCP , double *dparamLCP);

  /** lcp_cpg is a CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.
      \param nn      On enter, an integer which represents the dimension of the system.
      \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
      \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
      \param z       On return, a nn-vector of doubles which contains the solution of the problem.
      \param w       On return, a nn-vector of doubles which contains the solution of the problem.
      \param info    On return, an integer which contains the termination value:
      0: convergence\n
      1: iter = itermax\n
      2: negative diagonal term\n
      3: pWp nul

      \param iparamLCP On enter/return, a vector of integers:
      - iparamLCP[0] = itermax On enter, an integer which represents the maximum number of iterations allowed.
      - iparamLCP[1] = ispeak  On enter, an integer which represents the output log identifiant:
      0 : no output\n
      > 0: active screen ouput
      - iparamLCP[2] = it_end  On return, a double which represents the number of iterations performed by the algorithm.

      \param dparamLCP On enter/return, a vector of doubles:
      - dparamLCP[0] = tol     On enter, a double which represents the tolerance required.
      - dparamLCP[1] = res     On return, a double which represents the final error value.

      \author Mathieu Renouf.
  */
  void lcp_cpg(int *nn , double *M , double *q , double *z , double *w , int *info ,
               int *iparamLCP , double *dparamLCP);


  /** lcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for LCP.\n
      \param nn      On enter, an integer which represents the dimension of the system.
      \param M     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
      \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
      \param z       On return, a nn-vector of doubles which contains the solution of the problem.
      \param w       On return, a nn-vector of doubles which contains the solution of the problem.
      \param info    On return, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : negative diagonal term
      \param iparamLCP  On enter/return a vector of integers:\n
      - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
      - iparamLCP[1] = verbose  On enter, the output log identifiant:\n
      0 : no output\n
      >0: active screen output\n
      - iparamLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.
      \param dparamLCP  On enter/return a vector of doubles:\n
      - dparamLCP[0] = tol     On enter, the tolerance required.
      - dparamLCP[1] = omega   On enter, the relaxation parameter (not yet available).
      - dparamLCP[2] = res     On return, the final error value.
      \author Mathieu Renouf
  */
  void lcp_pgs(int *nn , double *M , double *q , double *z , double *w , int *info ,
               int *iparamLCP , double *dparamLCP);


  /** lcp_rpgs (Regularized Projected Gauss-Seidel ) is a solver for LCP, able to handle matrices with null diagonal terms.\n

  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
  0 : convergence\n
  1 : iter = itermax\n
  2 : negative diagonal term

  \param iparamLCP  On enter/return a vector of integers:\n
  - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
  - iparamLCP[1] = ispeak  On enter, the output log identifiant:\n
  0 : no output\n
  >0: active screen output\n
  - iparamLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.

  \param dparamLCP  On enter/return a vector of doubles:\n
  - dparamLCP[0] = tol     On enter, the tolerance required.
  - dparamLCP[1] = rho     On enter, the suggested regularization parameter
  - dparamLCP[2] = omega   On enter, the relaxation parameter (not yet available).
  - dparamLCP[3] = res     On return, the final error value.

  \author Mathieu Renouf & Pascal Denoyelle
  \todo Sizing the regularization paramter and apply it only on null diagnal term

  */
  void lcp_rpgs(int *nn , double *M , double *q , double *z , double *w , int *info ,
                int *iparamLCP , double *dparamLCP);

  /** lcp_psor Projected Succesive over relaxation solver for LCP. See cottle, Pang Stone Chap 5

  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
  0 : convergence\n
  1 : iter = itermax\n
  2 : negative diagonal term

  \param iparamLCP  On enter/return a vector of integers:\n
  - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
  - iparamLCP[1] = verbose  On enter, the output log identifiant:\n
  0 : no output\n
  >0: active screen output\n
  - iparamLCP[2] = initmethod On enter, the method of initialization
  0 : default w = q, z unchang

  it_end  On enter, the number of iterations performed by the algorithm.

  \param dparamLCP  On enter/return a vector of doubles:\n
  - dparamLCP[0] = tol     On enter, the tolerance required.
  - dparamLCP[1] = omega   On enter, the relaxation parameter
  - dparamLCP[2] = res     On return, the final error value.

  \author  Vincent Acary

  \todo use the relax parameter
  \todo add test
  \todo add a vector of relaxation parameter wtith an auto sizing (see SOR algorithm for linear solver.)

  */
  void lcp_psor(int *nn , double *M , double *q , double *z , double *w , int *info ,
                int *iparamLCP , double *dparamLCP);

  /** lcp_nsqp use a quadratic programm formulation for solving an non symmetric LCP
   *
   * Generic lcp parameters:\n
   *
   * \param nn      On enter, an integer which represents the dimension of the system.
   * \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
   * \param qq      On enter, a nn-vector of doubles which contains the components of the right hand side vector.
   * \param z       On return, a nn-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param w       On return, a nn-vector of doubles which contains the solution of the problem.
   * \param info    On return, an integer which returns the termination value:\n
   *                0 : convergence  / minimization sucessfull\n
   *                1 : Too Many iterations\n
   *            2 : Accuracy insuficient to satisfy convergence criterion\n
   *                5 : Length of working array insufficient\n
   *                Other : The constraints are inconstent\n
   *
   * Specific NSQP parameters:\n
   *
   * \param iparamLCP    On enter/return, a vector of integers (not used here).
   * \param dparamLCP    On enter/return, a vector of doubles:\n
   *                     - dparamLCP[0] = tol   On enter, the tolerance required.
   *
   * \author Vincent Acary
   *
   */
  void lcp_nsqp(int *nn , double *M , double *q , double *z , double *w , int *info ,
                int *iparamLCP , double *dparamLCP);

  /** lcp_latin (LArge Time INcrements) is a basic latin solver for LCP.

  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
  0 : convergence\n
  1 : iter = itermax\n
  2 : Cholesky Factorization failed \n
  3 : nul diagonal term\n

  \param iparamLCP On enter/return a vector of integers,
  - iparamLCP[0] = itermax  On enter, the maximum number of iterations allowed.
  - iparamLCP[1] = iout     On enter, the output log identifiant:\n
  0: no output\n
  >0: active screen output\n
  - iparamLCP[2] = it_end   On return, the number of iterations performed by the algorithm.
  \param dparamLCP  On enter/return a vector of doubles,
  - dparamLCP[0] = tol      On enter, the tolerance required.
  - dparamLCP[1] = k_latin  On enter, the latin parameter (a double strictly positive).
  - dparamLCP[2] = res      On return, the final error value.

  \author Nineb Sheherazade.
  */
  void lcp_latin(int *nn , double *M , double *q , double *z , double *w , int *info ,
                 int *iparamLCP , double *dparamLCP);

  /** lcp_latin_w (LArge Time INcrements) is a basic latin solver with relaxation for LCP.\n


  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
  0 : convergence\n
  1 : iter = itermax\n
  2 : Cholesky Factorization failed \n
  3 : nul diagonal term\n

  \param iparamLCP  On enter/return, a vector of integers:
  - iparamLCP[0] = itermax  On enter, the maximum number of iterations allowed.
  - iparamLCP[1] = chat     On enter, the output log identifiant:\n
  0 : no output\n
  >0 : active screen output\n
  - iparamLCP[2] = it_end   On return, the number of iterations performed by the algorithm.

  \param dparamLCP  On enter/return, a vector of doubles:
  - dparamLCP[0] = tol    On enter, the tolerance required.
  - dparamLCP[1] = k_latin  On enter, the latin parameter( a double strictly positive).
  - dparamLCP[2] = res      On return, the final error value.
  - dparamLCP[3] = relax    On enter, the relax coefficient ( \f$0 \leq relax \leq 1 \f$).

  \author Nineb Sheherazade.
  */
  void lcp_latin_w(int *nn , double *M , double *q , double *z , double *w , int *info ,
                   int *iparamLCP , double *dparamLCP);

  /** lcp_lexicolemke is a direct solver for LCP based on pivoting method principle for degenerate problem.\n
      Choice of pivot variable is performed via lexicographic ordering
      Ref: "The Linear Complementary Problem" Cottle, Pang, Stone (1992)\n


      \param nn      On enter, an integer which represents the dimension of the system.
      \param vec     On enter, a (\f$nn\times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
      \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
      \param zlem    On return, a nn-vector of doubles which contains the solution of the problem.
      \param wlem    On return, a nn-vector of doubles which contains the solution of the problem.
      \param info    On return, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : negative diagonal term\n

      \param iparamLCP  On enter/return, a vetor of integers:\n
      - iparamLCP[0] = itermax On enter, the maximum number of pivots allowed.
      - iparamLCP[1] = ispeak  On enter, the output log identifiant:\n
      0 : no output\n
      >0: active screen output\n
      - iparamLCP[2] = it_end  On return, the number of pivots performed by the algorithm.

      \param dparamLCP  On enter/return, a vetor of doubles (not used).\n



      \author Mathieu Renouf

  */
  void lcp_lexicolemke(int *nn , double *M , double *q , double *z , double *w , int *info ,
                       int *iparamLCP , double *dparamLCP);

  /** lcp_newton_min uses a nonsmooth Newton method based on the min formulation  (or max formulation) of the LCP

  \f$
  0 \le z \perp w \ge 0 \Longrightarrow \min(w,\rho z)=0 \Longrightarrow w = \max(0,w - \rho z)
  \f$

  \f$
  H(z) = H(\left[ \begin{array}{c} z \\ w \end{array}\right])= \left[ \begin{array}{c} w-Mz-q \\ min(w,\rho z) \end{array}\right] =0\\
  \f$


  References: Alart & Curnier 1990, Pang 1990


  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
  0 : convergence\n
  1 : iter = itermax\n
  2 : Problem in resolution in DGESV\n

  \param iparamLCP On enter/return a vectr of integers:\n
  - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
  - iparamLCP[1] = verbose  On enter, the output log identifiant:\n
  0 : no output\n
  >0 : active screen output\n
  - iparamLCP[2] = it_end  On return, the number of iterations performed by the algorithm.

  \param dparamLCP On enter/return, a vector of doubles:\n
  - dparamLCP[0] = tol     On enter the tolerance required.
  - dparamLCP[1] = res     On return the final error value.

  \author Vincent Acary

  \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
  \todo Add rules for the computation of the penalization rho
  \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
  */
  void lcp_newton_min(int *nn , double *M , double *q , double *z , double *w , int *info ,
                      int *iparamLCP , double *dparamLCP);


  /** lcp_newton_min use a nonsmooth newton method based on the Fischer-Bursmeister convex function
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
   * Generic lcp parameters:\n
   *
   * \param nn      Unchanged parameter which represents the dimension of the system.
   * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param info    Modified parameter which returns the termination value\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *                2 - negative diagonal term\n
   *
   * Specific NLGS parameters:\n
   *
   * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
   * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
   *                       0 - no output\n
   *                       0 < active screen output\n
   * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
   *
   * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
   * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
   *
   * \author Vincent Acary
   *
   * \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
   * \todo Add rules for the computation of the penalization rho
   * \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
   */
  void lcp_newton_FB(int *nn , double *M , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP);

  /**  lcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for LCP.
       \param nn      On enter, an integer which represents the dimension of the system.
       \param M     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
       \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
       \param z       On return, a nn-vector of doubles which contains the solution of the problem.
       \param w       On return, a nn-vector of doubles which contains the solution of the problem.
       \param info    On return, an integer which returns the termination value:\n
       0 : convergence\n
       1 : iter = itermax\n
       2 : negative diagonal term
       \param iparamLCP  On enter/return a vector of integers:\n
       - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
       - iparamLCP[1] = verbose  On enter, the output log identifiant:\n
       0 : no output\n
       >0: active screen output\n
       - iparamLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.
       \param dparamLCP  On enter/return a vector of doubles:\n
       - dparamLCP[0] = tol     On enter, the tolerance required.
       - dparamLCP[1] = omega   On enter, the relaxation parameter (not yet available).
       - dparamLCP[2] = res     On return, the final error value.
       \author Olivier Bonnefon
  */
  void lcp_path(int *nn , double *M , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

  /**
  * This function checks the validity of the vector z as a solution \n
  * of the LCP : \n
  * \f$
  *    0 \le z \perp Mz + q \ge 0
  * \f$
  * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
  * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
  * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
  * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
  * \author Pascal Denoyelle
  */
  int filter_result_LCP(int n, double *M , double *q , double *z , double tol, int chat, double *w);

  /**
   * This function checks the validity of the vector z as a solution \n
   * of the LCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
   * \author Vincent Acary form the routine  filter_result_LCP.c of Pascal Denoyelle
   */
  int lcp_compute_error(int n, double *M , double *q , double *z , int chat, double *w, double * error);

  /** Function to check convergence after LCP computation
  * \author Franck Perignon
   */
  int filter_result_LCP_block(SparseBlockStructuredMatrix *blmat, double *q , double *z , double tol, int chat, double *w);

#ifdef __cplusplus
}
#endif

#endif

