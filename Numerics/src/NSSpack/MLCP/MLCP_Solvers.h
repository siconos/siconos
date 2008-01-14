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
#ifndef MLCP_SOLVERS_H
#define MLCP_SOLVERS_H

/*! \page MLCPSolvers Mixed Linear Complementary Problems (MCLP)
  \section mlcpIntro The problem
  Try \f$(u,v,w)\f$ such that:\n
  \f$
  \left\lbrace
  \begin{array}{l}
  A u + Cv +a =0\\
  D u + Bv +b = w \\
  0 \le v \perp  w \ge 0\\
  \end{array}
  \right.
  \f$

  where  A is an (\f$ n \times n\f$ ) matrix, B is an (\f$ m \times m\f$ ) matrix,  C is an (\f$ n \times m\f$ ) matrix,\n
  D is an (\f$ m \times n\f$ ) matrix,    a and u is an (\f$ n \f$ ) vectors b,v and w is an (\f$ m \f$ ) vectors.

  \section mlcpSolversList Available solvers
  Use the generic function mlcp_solver(), to call one the the specific solvers listed below:
  - mlcp_pgs(), (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP
  - mlcp_rpgs(), (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP
  - mlcp_psor(), (projected successive overrelaxation method) is a solver for MLCP
  - mlcp_rpsor(),(regularized projected successive overrelaxation method) is a solver for MLCP
  - mlcp_path(), path solver

  The structure method, argument of mlcp_solver(), is used to give the name and parameters of the required solver.

  (see the functions/solvers list in mlcp_solvers.h)

*/

/*!\file MLCP_Solvers.h
  Solvers for Mixed Linear Complementary Problems (MCLP)
  \author Vincent Acary
*/


/** A type definition for a structure method_mlcp.
    \param name       name of the solver.
    \param itermax    maximum number of iterations.
    \param tol        convergence criteria value.
    \param relax      relaxation coefficient.
    \param rho        regularization coefficient
    \param chat       output boolean ( 0 = no output log ).
    \param iter       final number of iterations.
    \param err        final value of error criteria.
*/
typedef struct
{
  char   name[64];
  int    itermax;
  double tol;
  double relax;
  double rho;
  int    chat;
  int    iter;
  double err;
} method_mlcp;

#ifdef __cplusplus
extern "C" {
#endif

  /** mlcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
      \param A            On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
      \param B            On enter, a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
      \param C            On enter, a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
      \param D            On enter, a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
      \param a            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
      \param b            On enter, a m-vector of doubles which contains the components of the constant right hand side vector.
      \param nn           On enter, an integer which represents one the dimension of the MLCP problem.
      \param mm           On enter, an integer which represents one the dimension of the MLCP problem.
      \param u            On return, a n-vector of doubles which contains the solution of the problem.
      \param v            On return, a m-vector of doubles which contains the solution of the problem.
      \param w            On return, a m-vector of doubles which contains the complementary solution of the problem.
      \param info    On return, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : negative diagonal term
      \param iparamMLCP  On enter/return a vector of integers:\n
      - iparamMLCP[0] = itermax On enter, the maximum number of iterations allowed.
      - iparamMLCP[1] = verbose  On enter, the output log identifiant:\n
      0 : no output\n
      >0: active screen output\n
      - iparamMLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.
      \param dparamMLCP  On enter/return a vector of doubles:\n
      - dparamMLCP[0] = tol     On enter, the tolerance required.
      - dparamMLCP[1] = omega   On enter, the relaxation parameter (not yet available).
      - dparamMLCP[2] = res     On return, the final error value.
      \author Vincent Acary
  */
  void mlcp_pgs(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  /** mlcp_rpgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.
      \param A            On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
      \param B            On enter, a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
      \param C            On enter, a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
      \param D            On enter, a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
      \param a            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
      \param b            On enter, a m-vector of doubles which contains the components of the constant right hand side vector.
      \param nn           On enter, an integer which represents one the dimension of the MLCP problem.
      \param mm           On enter, an integer which represents one the dimension of the MLCP problem.
      \param u            On return, a n-vector of doubles which contains the solution of the problem.
      \param v            On return, a m-vector of doubles which contains the solution of the problem.
      \param w            On return, a m-vector of doubles which contains the complementary solution of the problem.
      \param info    On return, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : negative diagonal term
      \param iparamMLCP  On enter/return a vector of integers:\n
      - iparamMLCP[0] = itermax On enter, the maximum number of iterations allowed.
      - iparamMLCP[1] = verbose  On enter, the output log identifiant:\n
      0 : no output\n
      >0: active screen output\n
      - iparamMLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.
      \param dparamMLCP  On enter/return a vector of doubles:\n
      - dparamMLCP[0] = tol     On enter, the tolerance required.
      - dparamMLCP[1] = omega   On enter, the relaxation parameter (not yet available).
      - dparamMLCP[2] = res     On return, the final error value.
      \author Vincent Acary
  */
  void mlcp_rpgs(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  /** mlcp_psor (projected successive overrelaxation method) is a solver for MLCP.
      \param A            On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
      \param B            On enter, a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
      \param C            On enter, a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
      \param D            On enter, a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
      \param a            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
      \param b            On enter, a m-vector of doubles which contains the components of the constant right hand side vector.
      \param nn           On enter, an integer which represents one the dimension of the MLCP problem.
      \param mm           On enter, an integer which represents one the dimension of the MLCP problem.
      \param u            On return, a n-vector of doubles which contains the solution of the problem.
      \param v            On return, a m-vector of doubles which contains the solution of the problem.
      \param w            On return, a m-vector of doubles which contains the complementary solution of the problem.
      \param info    On return, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : negative diagonal term
      \param iparamMLCP  On enter/return a vector of integers:\n
      - iparamMLCP[0] = itermax On enter, the maximum number of iterations allowed.
      - iparamMLCP[1] = verbose  On enter, the output log identifiant:\n
      0 : no output\n
      >0: active screen output\n
      - iparamMLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.
      \param dparamMLCP  On enter/return a vector of doubles:\n
      - dparamMLCP[0] = tol     On enter, the tolerance required.
      - dparamMLCP[1] = omega   On enter, the relaxation parameter (not yet available).
      - dparamMLCP[2] = res     On return, the final error value.
      \author Vincent Acary
  */
  void mlcp_psor(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  /** mlcp_rpsor (regularized projected successive overrelaxation method) is a solver for MLCP.
      \param A            On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
      \param B            On enter, a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
      \param C            On enter, a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
      \param D            On enter, a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
      \param a            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
      \param b            On enter, a m-vector of doubles which contains the components of the constant right hand side vector.
      \param nn           On enter, an integer which represents one the dimension of the MLCP problem.
      \param mm           On enter, an integer which represents one the dimension of the MLCP problem.
      \param u            On return, a n-vector of doubles which contains the solution of the problem.
      \param v            On return, a m-vector of doubles which contains the solution of the problem.
      \param w            On return, a m-vector of doubles which contains the complementary solution of the problem.
      \param info    On return, an integer which returns the termination value:\n
      0 : convergence\n
      1 : iter = itermax\n
      2 : negative diagonal term
      \param iparamMLCP  On enter/return a vector of integers:\n
      - iparamMLCP[0] = itermax On enter, the maximum number of iterations allowed.
      - iparamMLCP[1] = verbose  On enter, the output log identifiant:\n
      0 : no output\n
      >0: active screen output\n
      - iparamMLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.
      \param dparamMLCP  On enter/return a vector of doubles:\n
      - dparamMLCP[0] = tol     On enter, the tolerance required.
      - dparamMLCP[1] = omega   On enter, the relaxation parameter (not yet available).
      - dparamMLCP[2] = res     On return, the final error value.
      \author Vincent Acary
  */
  void mlcp_rpsor(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  /** path solver
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
  void mlcp_path(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

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
  int mlcp_filter_result(int* n, int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v,  double tol, int chat, double *w);

  /**
   * This function checks the validity of the vector z as a solution \n
   * of the MLCP : \n
   * \f$
   *  \left\lbrace
   *   \begin{array}{l}
   *   A u + Cv +a =0\\
   *   D u + Bv +b = w
   *   0 \le v \perp  w \ge 0\\
   *   \end{array}
   *  \right.
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * It changes the input vector w by storing \f$ Mz + q \f$ in it.\n
   * \author Vincent Acary form the routine  filter_result_LCP.c of Pascal Denoyelle
   */
  int mlcp_compute_error(int* n, int* mm,  double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v,  int chat, double *w, double * error);

#ifdef __cplusplus
}
#endif

#endif
