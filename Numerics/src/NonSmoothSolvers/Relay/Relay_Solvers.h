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

/*! \page RelaySolvers Relay problems
  \section relayIntro
  Two formulations are available, primal and dual problems.

  Try \f$(z,w)\f$ such that:\n

  <em>primal problem:</em>\n

  \f$
  \left\lbrace
  \begin{array}{l}
  w - M z = q \\
  -w \in \partial\psi_{[-b, a]}(z)\\
  \end{array}
  \right.
  \f$

  or

  <em>dual problem:</em>\n
  \f$
  \left\lbrace
  \begin{array}{l}
  w - M z = q\\
  -z \in \partial\psi_{[-b,a]}(w)\\
  \end{array}
  \right.
  \f$

  here M is an (\f$ n \times n \f$)-matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

  Primal case:
  Use the generic function pr_solver(), to call one the the specific solvers listed below:

  - pr_latin(), latin solver for primal relay problems.
  - pr_nlgs(), non linear Gauss-Seidel solver for primal relay problems

  Dual case:
  Use the generic function dr_solver(), to call one the the specific solvers listed below:

  - dr_latin(), latin (LArge Time INcrement)solver for dual relay problems.
  - dr_nlgs(), non linear Gauss-Seidel solver for dual relay problems

  The structure method, argument of pr_solver() or dr_solver(), is used to give the name and parameters of the required solver.

  (see the functions/solvers list in Relay_solvers.h)

*/
/*!\file Relay_Solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
  Subroutines for the resolution of relay problems.
*/

/** A type definition for a structure method_pr.
    \param name       name of the solver.
    \param itermax    maximum number of iterations.
    \param tol        convergence criteria value.
    \param k_latin    latin coefficient
    \param a          upper bound.
    \param b          lower bound.
    \param chat       output boolean ( 0 = no output log ).
    \param normType   name norm (not yet available).

    \param iter       final number of iterations.
    \param err        final value of error criteria.
*/
typedef struct
{

  char     name[64];
  int      itermax;
  double   tol;
  double   k_latin;
  double   *a;
  double   *b;
  int      chat;
  char     normType[64];
  int      iter;
  double   err;

} method_pr;

/** A type definition for a structure method_dr.

\param name       name of the solver.
\param itermax    maximum number of iterations.
\param tol        convergence criteria value.
\param k_latin    latin coefficient
\param a          upper bound
\param b          lower bound
\param chat       output boolean ( 0 = no output log ).
\param normType   name norm (not yet available).

\param iter       final number of iterations
\param err        final value of error criteria
*/
typedef struct
{

  char     name[64];
  int      itermax;
  double   tol;
  double   k_latin;
  double   *a;
  double   *b;
  int      chat;
  char     normType[64];
  int      iter;
  double   err;

} method_dr;

#ifdef __cplusplus
extern "C" {
#endif

  /** pr_latin is a specific latin solver for primal relay problems.


  \param vec       On enter, a (\f$ n \times n\f$)-vector of doubles containing the components of the matrix with a fortran storage.
  \param qq        On enter, a n-vector of doubles containing the components of the vector.
  \param n        On enter, an integer which represents the dimension of the second member.
  \param k_latin   On enter, a double, the latin coefficient (strictly non negative).
  \param a         On enter, a n-vector of doubles, the upper bound.
  \param b         On enter, a n-vector of doubles, the down bound.
  \param itermax   On enter, an integer which represents the maximum iterations required.
  \param tol       On enter, a double which contains the tolerance required.
  \param chat      On enter, an integer the output log identifiant:\n                              0 : no output\n
  >0 : active screen output
  \param it_end    On return, an integer which represents the number of iterations carried out.
  \param res       On return, a double, the error value.
  \param z         On return,a n-vector of doubles wich contains the solution of the problem.
  \param w         On return, a n-vector of doubles which contains the solution of the problem.
  \param info      On return, an integer which represents the termination reason: \n
  0 = convergence,\n
  1 = no convergence,\n
  2 = Cholesky factorization failed,\n
  3 = Nul diagonal term\n


  \author Nineb Sheherazade.
  */
  void pr_latin(double* , double* , int* , double* , double* , double* , int* , double* , int *, double* , double* , int* , double* , int*);

  /** pr_nlgs is a specific nlgs(non linear Gauss-Seidel) solver for primal relay problems.\n



  \param vec       On enter, a (\f$ n\times n \f$)-vector of doubles containing the components of the matrix with a fortran storage.
  \param qq        On enter, a n-vector of doubles containing the components of the vector.
  \param n        On enter, an integer, the dimension of the second member.
  \param a         On enter, a n-vector of doubles, the upper bound.
  \param b         On enter, a n-vector of doubles, the down bound.
  \param itermax   On enter, an integer, the maximum iterations required.
  \param tol       On enter, a double, the tolerance required.
  \param chat      On enter, an integer, the output log identifiant:\n
  0 =  no output \n
  >0 =  active screen output\n


  \param it_end    On return, an integer, the number of iterations carried out.
  \param res       On return, a double, the error value.
  \param z         On return, a n-vector of doubles, the solution of the problem.
  \param w         On return, a n-vector of double, the solution of the problem.
  \param info      On return, an integer, the termination reason:\n
  0 = convergence,\n
  1 = no convergence,\n
  2 = Nul diagonal term\n


  \author Nineb Sheherazade.
  */
  void pr_nlgs(double* , double* , int* , double* , double* , int* , double* , int*, double* , double* , int* , double* , int *);

  /** dr_latin is a specific latin (LArge Time INcrement)solver for dual relay problems.\n

  \param vec        On enter, a (\f$ n\times n \f$)-vector of doubles containing the components of the double matrix with a fortran storage.
  \param q          On enter, a n-vector of doubles containing the components of the second member.
  \param n         On enter, an integer, the dimension of the second member.
  \param a          On enter, a n-vector of doubles, the upper bound.
  \param b          On enter, a n-vector of doubles, the lower bound.
  \param itermax    On enter, an integer, the maximum iterations required.
  \param tol        On enter, a double, the tolerance required.
  \param k_latin    On enter, a double, the search direction (strictly non negative).
  \param chat       On enter, an integer, the output log identifiant:\n
  0  =  no output \n
  >0 =  active screen output\n

  \param it_end     On return, an integer, the number of iterations carried out.
  \param res        On return, a double, the error value.
  \param z          On return, a n-vector of doubles, the solution of the problem.
  \param w          On return, a n-vector of doubles, the solution of the problem.
  \param info       On return, an integer, the termination reason:\n
  0 = convergence,\n
  1 = non convergence,\n
  2 = Cholesky factorization failed,\n
  3 = Nul diagonal term.\n

  \author Nineb Sheherazade.
  */
  void dr_latin(double * , double *, int *, double * , double *, double *, int *, double *, int *, double* , double* , int *, double *, int *)  ;

  /**  dr_nlgs is a specific nlgs (Non Linear Gauss Seidel) solver for dual relay problems.\n

  \param vec      On enter, a (\f$ n\times n \f$)-vector of doubles containing the components of the double matrix with a fortran90 allocation.
  \param q        On enter, a n-vector of doubles containing the components of the double vector.
  \param n       On enter, an integer, the dimension of the second member.
  \param a        On enter, a n-vector of doubles, the upper bound.
  \param b        On enter, a n-vector of doubles, the lower bound.
  \param itermax  On enter, an integer, the maximum iterations required.
  \param tol      On enter, a double, the tolerance required.
  \param chat     On enter, an integer, the output log identifiant:\n
  0  =  no output, \n
  >0 =  active screen output\n

  \param it_end   On return, an integer, the number of iterations carried out.
  \param res      On return, a double, the error value.
  \param z        On return, a n-vector of doubles, the solution of the problem.
  \param w        On return, a n-vector of doubles, the solution of the problem.
  \param info     On return, an integer, the termination reason:\n
  0 = convergence,\n
  1 = no convergence,\n
  2 = nul diagonal term.\n


  \author Nineb Sheherazade.
  */
  void dr_nlgs(double *vec , double *q , int *n , double *a , double *b , int *itermax , double *tol , int *chat,
               double *z , double *w , int *it_end , double *res , int *info);

  /** pr_gsnl is a specific gsnl (Gauss Seidel Non Linear)solver for relay problems.
      \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
      \param qq On enter a pointer over doubles containing the components of the double vector.
      \param nn On enter a pointer over integers, the dimension of the second member.
      \param a On enter a pointer over doubles, the bound.
      \param itermax On enter a pointer over integers, the maximum iterations required.
      \param tol On enter a pointer over doubles, the tolerance required.
      \param it_end On enter a pointer over integers, the number of iterations carried out.
      \param res On return a pointer over doubles, the error value.
      \param z On return double vector, the solution of the problem.
      \param w On return double vector, the solution of the problem.
      \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
      \author Nineb Sheherazade.
  */
  void pr_gsnl(double vec[], double *q, int *nn, double a[], double b[], int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info);

#ifdef __cplusplus
}
#endif

#endif
