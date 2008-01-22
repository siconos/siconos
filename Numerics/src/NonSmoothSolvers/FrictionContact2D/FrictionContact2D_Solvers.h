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
#ifndef FrictionContact2DSolvers_H
#define FrictionContact2DSolvers_H

/*! \page fc2DSolvers Friction-contact problems (2-dimensional)
  Resolution of contact problems with friction in the 2 dimensional case
  \section fc3DIntro The problem
  Two formulations are available, primal and dual problems.

  Try \f$(reaction,velocity)\f$ such that:

  <em>primal problem:</em>\n
  \f$
  \left\lbrace
  \begin{array}{l}
  velocity - M reaction = q \\
  0 \le reaction_n \perp velocity_n \ge 0\\
  -velocity_t \in \partial\psi_{[-\mu reaction_n, \mu reaction_n]}(reaction_t)\\
  \end{array}
  \right.
  \f$

  or \n

  <em>dual problem:</em>\n
  \f$
  \left\lbrace
  \begin{array}{l}
  velocity - M reaction = q \\
  0 \le reaction_n \perp velocity_n \ge 0 \\
  -reaction_t \in \partial\psi_{[-\mu velocity_n, \mu velocity_n]}(velocity_t)\\
  \end{array}
  \right.
  \f$

  M is an (\f$ n\times n \f$)-matrix, q, reaction and velocity some n-dimensional vectors.

  \section fc2DSolversList Available solvers
  <em> For primal problem: </em>\n
  Use the generic function pfc_2D_solver(), to call one the the specific solvers listed below:
  - pfc_2D_latin(), latin solver
  - pfc_2D_nlgs(), Non Linear Gauss Seidel solver
  - pfc_2D_cpg(), conjugated projected gradient solver

  <em>For dual problem: </em>\n
  Use the generic function dfc_2D_solver(), to call one the the specific solvers listed below:
  - dfc_2D_latin(), latin solver

  The structure method, argument of pfc_2D_solver() or dfc_2D_solver(), is used to give the name and parameters of the required solver.
  (see the functions/solvers list in FrictionContact2D_solvers.h)

*/

/*!\file FrictionContact2D_Solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
  Subroutines for the resolution of contact problems with friction (2-dimensional case).\n
*/


/** A type definition for a structure method_pfc_2D.
    \param name       name of the solver.
    \param itermax    maximum number of iterations.
    \param tol        convergence criteria value.
    \param k_latin    search direction of the latin metod.
    \param chat       output boolean ( 0 = no output log ).
    \param normType   name norm (not yet available).
    \param iter       final number of iterations.
    \param err        final value of error criteria
*/
typedef struct
{
  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  int    chat;
  char   normType[64];
  int    iter;
  double err;

} method_pfc_2D;

/** A type definition for a structure method_dfc_2D
    \param name       name of the solver.
    \param itermax    maximum number of iterations.
    \param normType   name norm (not yet available).
    \param tol        convergence criteria value.
    \param k_latin    latin coefficient
    \param J1         gap in normal contact direction.
    \param ddl_n      the contact in normal direction dof (not prescribed),
    \param ddl_tt     the contact in tangential direction dof (not prescribed)
    \param ddl_d      the prescribed dof.
    \param dim_tt     the dimension of the vector ddl_tt.
    \param dim_d      the dimension of the vector ddl_d.
    \param chat       output boolean ( 0 = no output log ).
    \param iter       final number of iteration
    \param err        final value of error criteria
*/
typedef struct
{

  char   name[64];
  int    itermax;
  char   normType[64];
  double tol;
  double k_latin;
  double *J1;
  int    *ddl_n;
  int    *ddl_tt;
  int    *ddl_d;
  int    dim_tt;
  int    dim_d;

  int   chat;
  int    iter;
  double err;

} method_dfc_2D;
#ifdef __cplusplus
extern "C" {
#endif

  /** specific cpg (conjugated projected gradient) solver for primal contact problems with friction in the 2D case.
      \param nc          On enter an integer, the number of contacst. The dimension of the system is n=2*nc.
      \param vec         On enter a (\f$ n\times n \f$)-vector of doubles containing the components of the double matrix with a fortran90 allocation ( M ).
      \param b           On enter a n-vector of doubles containing the components of the double vector ( q ).
      \param mu          On enter, the vector of friction coefficients (mu[i] corresponds to contact i)
      \param iparamPFC   On enter/return a vector of integers:\n
          - iparamPFC[0] = on enter, the maximum number of iterations allowed,
          - iparamPFC[1] = on enter, the parameter which represents the output log identifiant,\n
             0 - no output\n
             >0 - active screen output\n
          - iparamPFC[2] =  on return, the number of iterations performed by the algorithm.
      \param dparamPFC    On enter/return a vector of doubles:\n
          - dparamPFC[0] = on enter, a positive double which represents the tolerance required,
          - dparamPFC[1] = on return, a positive double which represents the residu.
      \param x           On return a n-vector of doubles, the solution of the problem ( z ).
      \param rout        On return a n-vector of doubles, the solution of the problem ( w ).
      \param info        On return an integer, the termination reason:\n
          0 = convergence,\n
          1 = no convergence,\n
          2 = Operation of alpha no conform.\n
      \author Nineb Sheherazade.
  */
  void pfc_2D_cpg(int, double *vec , double *q , double *reaction , double *velocity , double *mu , int *info , int *iparamLCP , double *dparamLCP);

  /**   pfc_2D_gsnl is a specific nlgs (Non Linear Gauss Seidel ) solver for primal contact problem with friction in 2D case.
  \param nc          On enter an integer, the number of contacst. The dimension of the system is n=2*nc.
  \param vec         On enter a (\f$ n\times n \f$)-vector of doubles containing the components of the double matrix with a fortran90 allocation.
  \param qq          On enter a nn-vector of doubles containing the components of the second member.
  \param mu          On enter, the vector of friction coefficients (mu[i] corresponds to contact i)
  \param iparamPFC   On enter/return a vector of integers:\n
  - iparamPFC[0] = on enter, the maximum number of iterations allowed,
  - iparamPFC[1] = on enter, the parameter which represents the output log identifiant:\n
  0 - no output\n
  >0 - active screen output\n
  - iparamPFC[2] =  on return, the number of iterations performed by the algorithm.\n
  \param dparamPFC    On enter/return a vector of doubles:\n
  - dparamPFC[0] = on enter, a positive double which represents the tolerance required,
  - dparamPFC[1] = on return, a positive double which represents the residu.
  \param z           On return a nn-vector of doubles, the solution of the problem.
  \param w           On return a nn-vector of doubles, the solution of the problem.
  \param info        On return an integer, the termination reason:\n
  0 = convergence,\n
  1 = no convergence,\n
  2 = nul term in denominator.
  \author Nineb Sheherazade.
  */
  void pfc_2D_nlgs(int, double *vec , double *q , double *reaction , double *velocity, double *mu  , int *info , int *iparamLCP , double *dparamLCP);

  /**   pfc_2D_latin  is a specific latin solver for primal contact problem with friction in the 2D case.
  \param nc          On enter an integer, the number of contacst. The dimension of the system is n=2*nc.
  \param vec         On enter a ( \f$ n\times n \f$)-vector of doubles containing the components of the double matrix with a fortran90 allocation.
  \param qq          On enter a n-vector of doubles containing the components of the second member.
  \param mu          On enter, the vector of friction coefficients (mu[i] corresponds to contact i)
  \param iparamPFC   On enter/return a vector of integers:\n
  - iparamPFC[0] = on enter, the maximum number of iterations allowed,\n
  - iparamPFC[1] = on enter, the parameter which represents the output log identifiant:\n
  0 - no output\n
  >0 -  active screen output\n
  - iparamPFC[2] =  on return, the number of iterations performed by the algorithm.\n
  \param dparamPFC   On enter/return a vector of doubles:\n
  - dparamPFC[0] = on enter, a positive double which represents the tolerance required,
  - dparamPFC[1] = on enter, a strictly nonnegative double which represents the search parameter,\n
  - dparamPFC[2] = on return, a positive double which represents the residu.
  \param z           On return a n-vector of doubles, the solution of the problem.
  \param w           On return a n-vector of doubles, the solution of the problem.
  \param info        On return an integer, the termination reason:
  0 = Convergence,\n
  1 = no convergence,\n
  2 = Cholesky factorizayion failed,\n
  3 = Nul term in diagonal of M.\n
  \author Nineb Sheherazade.
  */
  void pfc_2D_latin(int, double *, double *, double *, double *, double *, int *, int *, double *);

  /** pfc_2D_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm for primal contact problem with friction.\n
   *
   *
   * \param xi        On enter, the intermediate iterate which goes to be projected (projc1).
   * \param nn        On enter, the dimension of the system.
   * \param statusi   On enter, a vector which contains the initial status.
   * \param pi        On enter, a vector which contains the components of the descent direction.
   * \param fric      On enter, a vector which contains the friction coefficient.
   * \param projc1    On return, the corrected iterate.
   * \param projc2    On return, the new status.
   *
   * \author Sheherazade Nineb.
   *
   */
  void pfc_2D_projc(int nc , double mu , double *reaction , double *p , int *status);

  /** pfc_2D_projf is a specific projection operator related to CPG (conjugated projected gradient) algorithm
   *              for primal contact problem with friction.\n
   *
   *
   * \param etat       On enter,  parameter which represents the status vector.
   * \param nn         On enter,  parameter which represents the dimension of the system.
   * \param y          On enter,  parameter which contains the components of the residue or descent direction vector.
   * \param fric       On enter,  parameter which contains the friction coefficient.
   * \param projf1     On return, parameter which contains the projected residue or descent direction.
   *
   * \author Shéhérazade Nineb.
   *
   */
  void pfc_2D_projf(int n , double *ww , double *zz , double *rr , double *pp , int *status);


  /**   cfd_latin  is a specific latin solver for dual contact problem with friction in the 2D case.\n
       \param vec      On enter a (\f$ n \times n\f$)-vector of doubles containing the components of the double matrix with a fortran storage.
       \param qq       On enter a nn-vector of doubles containing the components of the second member.
       \param nn       On enter an integer, the dimension of the second member.
       \param k_latin  On enter a double, the latin coefficient (strictly nonnegative).
       \param mu     On enter a positive double, the friction coefficient.
       \param itermax  On enter an integer, the maximum iterations required.
       \param tol      On enter a double, the tolerance required.
       \param z        On return a nn-vector of doubles, the solution of the problem.
       \param w        On return a nn-vector of doubles, the solution of the problem.
       \param it_end   On enter an integer, the number of iterations carried out.
       \param res      On return a double, the error value.
       \param info     On return an integer, the termination reason:\n
       0 = successful, \n
       1 = no convergence, \n
       2 = Cholesky failed, \n
       3 = nul diagonal term, \n
       4 = nul diagonal term, \n
       \author Nineb Sheherazade.
  */
  void dfc_2D_latin(double* , double* , int* , double* , double* , int* , double* , int *, double* , double* , int* , double* , int*);

  /**   This subroutine allows the formulation in the LCP (Linear  Complementary Problem) form
       of a 2D contact problem with friction.\n
       \param dim_F1    On enter a pointer over integers, the dimension of the DFC_2D problem,
       \param mu      On enter a pointer over doubles, the friction coefficient,
       \param K1        On enter a pointer over doubles containing the components of the
       rigidity matrix with a fortran90 storage,
       \param F1        On enter a pointer over doubles containing the right hand side,
       \param ddl_n     On enter a pointer over integers , the contact in normal direction dof
       (not prescribed),
       \param ddl_tt    On enter a pointer over integers, the contact in tangential direction dof
       (not prescribed)
       \param dim_nc    On enter a pointer over integers, the dimension of the vector ddl_tt.
       \param ddl_d     On enter a pointer over integers, the prescribed dof,
       \param dim_d     On enter a pointer over integers, the dimension of the vector ddl_d,
       \param J1        On enter a pointer over doubles, gap in normal contact direction.
       \n\n
       \param MM        On return a pointer over doubles containing the components of a double
       matrix (3*dim_nc,3*dim_nc) with a fortran90 allocation.
       \param q         On return a pointer over doubles, a double vector (3*dim_nc).
       \author Nineb Sheherazade.
  */
  void dfc_2D2lcp(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);

  /**   This routine allows to give the solution of the 2D contact problem with friction given. \n
       \param dim_F1    On enter a pointer over integers, the dimension of the DFC_2D problem,
       \param ztel      On enter a pointer over doubles, the solution given by a LCP solver.
       \param wtel      On enter a pointer over doubles, the solution given by a LCP solver.
       \param K1        On enter a pointer over doubles containing the components of the
       rigidity matrix with a fortran90 storage,
       \param F1        On enter a pointer over doubles containing the right hand side,
       \param J1        On enter a pointer over doubles, gap in normal contact direction.
       \param ddl_n     On enter a pointer over integers , the contact in normal direction dof
       (not prescribed),
       \param ddl_tt    On enter a pointer over integers, the contact in tangential direction dof
       (not prescribed)
       \param dim_tt    On enter a pointer over integers, the dimension of the vector ddl_tt.
       \param ddl_d     On enter a pointer over integers, the prescribed dof,
       \param dim_d     On enter a pointer over integers, the dimension of the vector ddl_d,
       \n\n
       \param U2        On return a pointer over doubles, the solution of the contact friction problem U2(dim_F1).
       \param F2        On return a pointer over doubles, the solution of the contact friction problem F2(dim_F1).
       \author Nineb Sheherazade.
  */
  void lcp2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                  int *, int *, int *,  double *, double *);

  /** Formulation in the condensed form of a 2D contact problem with friction (DFC_2D).\n
      \param dim_F1    On enter a pointer over integers, the dimension of the DFC_2D problem,
      \param mu      On enter a pointer over doubles, the friction coefficient,
      \param K1        On enter a pointer over doubles containing the components of the
      rigidity matrix with a fortran90 storage,
      \param F1        On enter a pointer over doubles containing the right hand side,
      \param ddl_n     On enter a pointer over integers , the contact in normal direction dof
      (not prescribed),
      \param ddl_tt    On enter a pointer over integers, the contact in tangential direction dof
      (not prescribed)
      \param dim_nc    On enter a pointer over integers, the dimension of the vector ddl_tt.
      \param ddl_d     On enter a pointer over integers, the prescribed dof,
      \param dim_d     On enter a pointer over integers, the dimension of the vector ddl_d,
      \param J1        On enter a pointer over doubles, gap in normal contact direction.
      \n\n
      \param MM        On return a pointer over doubles containing the components of a double
      matrix (2*dim_nc,2*dim_nc) with a fortran90 allocation.
      \param q         On return a pointer over doubles, a double vector (2*dim_nc).
      \author Nineb Sheherazade.
  */
  void dfc_2D2cond_2D(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);

  /** Solution of the 2D contact problem with friction given.\n
      \param dim_F1    On enter a pointer over integers, the dimension of the DFC_2D problem,
      \param ztel      On enter a pointer over doubles, the solution given by a dfc_2D solver.
      \param wtel      On enter a pointer over doubles, the solution given by a dfc_2D solver.
      \param K1        On enter a pointer over doubles containing the components of the
      rigidity matrix with a fortran90 storage,
      \param F1        On enter a pointer over doubles containing the right hand side,
      \param J1        On enter a pointer over doubles, gap in normal contact direction.
      \param ddl_n     On enter a pointer over integers , the contact in normal direction dof
      (not prescribed),
      \param ddl_tt    On enter a pointer over integers, the contact in tangential direction dof
      (not prescribed)
      \param dim_tt    On enter a pointer over integers, the dimension of the vector ddl_tt.
      \param ddl_d     On enter a pointer over integers, the prescribed dof,
      \param dim_d     On enter a pointer over integers, the dimension of the vector ddl_d,
      \n\n
      \param U2        On return a pointer over doubles, the solution of the contact friction problem U2(dim_F1).
      \param F2        On return a pointer over doubles, the solution of the contact friction problem F2(dim_F1).
      \author Nineb Sheherazade.
  */
  void cond_2D2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                      int *, int *, int *,  double *, double *);

  /** */
  void projf(int*, int*, double*, double*, double*);

  /** */
  void projc(double*, int*, int*, double*, double*, double*, int*);

#ifdef __cplusplus
}
#endif

#endif
